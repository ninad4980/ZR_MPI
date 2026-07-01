#include "Movers.h"

#include "PotentialSolver.h"



//====================================================
// PIC mover
//====================================================

void PICMover::move(double rfFactor)
{
    double lc[3];

    for (Particle &part : sp.particles)
    {
        if (part.STATUS != PARTICLE_ACTIVE)
            continue;

        //--------------------------------------------
        // logical coordinate
        //--------------------------------------------

        domain.Posv3d2lc(
            lc,
            part.pos
        );

        //--------------------------------------------
        // gather electric field
        //--------------------------------------------
//		vectorCyl ef_cyl;
//		ef_cyl.r=0.0; ef_cyl.theta=0.0; ef_cyl.zee=0.0;
//        ef_cyl =
//            domain.Estatic.gather(
//                lc,
//                0.0,
//                domain.dh[0]
//            );
		
		 vectorCyl ef_cyl =
		    domain.electricField(
    		    lc,
        		rfFactor
    			);
        //--------------------------------------------
        // gather magnetic field
        //--------------------------------------------
//		vectorCyl B_cyl;
//		B_cyl.r=0.0; B_cyl.theta=0.0; B_cyl.zee=0.0;
//        B_cyl =
//            domain.Bstatic.gather(
//                lc,
//                0.0,
//                domain.dh[0]
//            );
		
		vectorCyl B_cyl =
    		domain.magneticField(
        		lc,
        		rfFactor
    			);
        //--------------------------------------------
        // particle cylindrical position
        //--------------------------------------------

        vectorCyl pos_c =
            Cart2Cyl(part.pos);

        //--------------------------------------------
        // convert to Cartesian
        //--------------------------------------------

        vector3d E =
            cylV2cart(
                ef_cyl,
                pos_c
            );

        vector3d B =
            cylV2cart(
                B_cyl,
                pos_c
            );
        //--------------------------------------------
        // Boris push
        //--------------------------------------------
        phasespace ps =
            BorisPush(
                sp.dt,
                part.pos,
                part.vel,
                E,
                B,
                sp.charge,
                sp.mass
            );

        //--------------------------------------------
        // update particle
        //--------------------------------------------

        part.pos = ps.pos;

        part.vel = ps.vel;
    }

  

    //--------------------------------------------
    // MPI transfer
    //--------------------------------------------

    sp.applyDomainConditions();

    sp.transferParticles();
      //--------------------------------------------
    // boundaries
    //--------------------------------------------

    sp.checkBounds();
}


//====================================================
// DSMC mover
//====================================================

void DSMCMover::move(double)
{
    //--------------------------------------------
    // free flight
    //--------------------------------------------

    for (Particle &part : sp.particles)
    {
        if (part.STATUS != PARTICLE_ACTIVE)
            continue;

        part.pos =
            part.pos
            +
            part.vel * sp.dt;
    }

    //--------------------------------------------
    // wall interactions
    //--------------------------------------------

    sp.checkBounds();

    //--------------------------------------------
    // DSMC collisions
    //--------------------------------------------

    performCollisions();

    //--------------------------------------------
    // MPI transfer
    //--------------------------------------------

    sp.applyDomainConditions();

    sp.transferParticles();
}


void DSMCMover::performCollisions()
{
    //--------------------------------------------
    // future:
    //
    // Bird DSMC
    // null-collision
    // VHS/VSS
    //
    //--------------------------------------------
}


//====================================================
// Fluid mover
//====================================================

void FluidMover::move(double rfFactor)
{
    //-----------------------------------
    // 1. compute ion velocity
    //-----------------------------------

    computeVelocity( rfFactor);

    //-----------------------------------
    // 2. compute fluxes
    //-----------------------------------

    computeFluxes();
    //-----------------------------------
    // 3. continuity update
    //-----------------------------------

    advanceDensity();

    //-----------------------------------
    // 4. MPI sync
    //-----------------------------------

    fs.density.updateBoundaries(
        domain.mpi_rank,
        domain.mpi_size,
        1
    );
}


void FluidMover::computeVelocity(double rfFactor)
{
    int Nr = domain.Nr;
    int Nz = domain.Nz_phys;

    for (int k=0; k<Nz; k++)
    {
        int kg = k + 1;

        for (int i=0; i<Nr; i++)
        {
            int obj =
                domain.object[i + k*Nr];

            //----------------------------------
            // plasma + extraction
            //----------------------------------

            if (
                obj == OBJ_PLASMA
                ||
                obj == OBJ_EXTRACTION
            )
            {
                vectorCyl E =
    				domain.electricField(
        			i,
        			kg,
			        rfFactor
    				);

				double Er = E.r;
				double Ez = E.zee;

                fs.ur(i,kg) =
                    fs.mobility * Er;

                fs.uz(i,kg) =
                    fs.mobility * Ez;
            }
            else
            {
                fs.ur(i,kg) = 0.0;
                fs.uz(i,kg) = 0.0;
            }
        }
    }

    fs.ur.updateBoundaries(
        domain.mpi_rank,
        domain.mpi_size,
        1
    );

    fs.uz.updateBoundaries(
        domain.mpi_rank,
        domain.mpi_size,
        1
    );
}

void FluidMover::computeFluxes()
{
	   int Nr = domain.Nr;
    int Nz = domain.Nz_phys;

    for (int k=0; k<Nz; k++)
    {
        int kg = k+1;

        for (int i=0; i<Nr; i++)
        {
            fs.flux_r(i,kg) =
                fs.density(i,kg)
                *
                fs.ur(i,kg);

            fs.flux_z(i,kg) =
                fs.density(i,kg)
                *
                fs.uz(i,kg);
        }
    }

    fs.flux_r.updateBoundaries(
        domain.mpi_rank,
        domain.mpi_size,
        1
    );

    fs.flux_z.updateBoundaries(
        domain.mpi_rank,
        domain.mpi_size,
        1
    );
}

void FluidMover:: advanceDensity()
{
	int Nr = domain.Nr;
    int Nz_phys = domain.Nz_phys;
    int Nz_local = domain.Nz_local;
	double dt = fs.dt;

    double dr = domain.dh[0];

    double dz = domain.dh[2];

    //----------------------------------------
    // copy old density
    //----------------------------------------

    for (int k=0; k<domain.Nz_local; k++)
    {
        for (int i=0; i<domain.Nr; i++)
        {
            fs.density_old(i,k)
                =
                fs.density(i,k);
        }
    }

    //----------------------------------------
    // continuity equation
    //----------------------------------------

    for (int k=1;
         k<domain.Nz_local-1;
         k++)
    {
        for (int i=1;
             i<domain.Nr-1;
             i++)
        {
            double r =
                i * dr;
			double rp = r + 0.5*dr;
			double rm = r - 0.5*dr;

			double Frp;
			double Frm;
            //----------------------------------
            // radial divergence
            //----------------------------------

           double urp = 0.5*(fs.ur(i,k)+fs.ur(i+1,k));

			double nface;

			if(urp>0.0)
			    nface = fs.density(i,k);
			else
			    nface = fs.density(i+1,k);

			 Frp = urp*nface;
			
			double urm = 0.5*(fs.ur(i,k)+fs.ur(i-1,k));

			if(urm>0.0)
			    nface = fs.density(i-1,k);
			else
			    nface = fs.density(i,k);

			 Frm = urm*nface;
			
			double div_r =
				(
			    rp*Frp
    			-
			    rm*Frm
				)
					/
				(r*dr);
            //----------------------------------
            // axial divergence
            //----------------------------------
				double uz = fs.uz(i,k);

				double Fp,Fm;

				if(uz>0)
					{
				    Fp = fs.flux_z(i,k);
				    Fm = fs.flux_z(i,k-1);
					}
				else
					{
				    Fp = fs.flux_z(i,k+1);
				    Fm = fs.flux_z(i,k);
					}

			double 	div_z = (Fp-Fm)/dz;
            //----------------------------------
            // update
            //----------------------------------

			fs.density(i,k)
				=
				fs.density_old(i,k)
				-
				dt*(div_r + div_z)
				+
				dt*fs.source(i,k);
            //----------------------------------
            // positivity limiter
            //----------------------------------

            if (fs.density(i,k) < 0.0)
                fs.density(i,k) = 0.0;
        }
    }

    fs.density.updateBoundaries(
        domain.mpi_rank,
        domain.mpi_size,
        1
    );
}
