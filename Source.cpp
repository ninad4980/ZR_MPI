#include "Source.h"

using namespace Const;




void Source::injectNeutralInlet(
    double pressure, double temperature,
    double inletRadius, double zInject,
    double dt,double massFlow
	)
{
    //--------------------------------------------------
    // only inject on rank 0
    //--------------------------------------------------
    if (domain.mpi_rank != 0)
        return;
    //--------------------------------------------------
    // gas properties
    //--------------------------------------------------

    double m = sp.mass;

    double area =
        PI * inletRadius * inletRadius;

    //--------------------------------------------------
    // number density
    //--------------------------------------------------

    double n0 =
        pressure / (kB * temperature);

    //--------------------------------------------------
    // thermal speed
    //--------------------------------------------------

    double vth =
        sqrt(
            8.0 * kB * temperature
            /
            (PI * m)
        );

    //--------------------------------------------------
    // molecular flux
    //--------------------------------------------------

    double Gamma =
        0.25 * n0 * vth;

    //--------------------------------------------------
    // real particles entering
    //--------------------------------------------------

    double num_real;

    //--------------------------------------------------
    // CASE 1:
    // known mass flow
    //--------------------------------------------------

    if (massFlow > 0.0)
    {
        num_real =
            (
                massFlow * dt
            )
            /
            m;
    }

    //--------------------------------------------------
    // CASE 2:
    // pressure-driven estimate
    //--------------------------------------------------

    else
    {
        num_real =
            Gamma * area * dt;
    }

    //--------------------------------------------------
    // simulation particles
    //--------------------------------------------------

    int num_sim =
        int(num_real / sp.spwt);

    double remain =
        num_real/sp.spwt - num_sim;

    if (rnd() < remain)
        num_sim++;

    //--------------------------------------------------
    // thermal sigma
    //--------------------------------------------------

    double sigma =
        sqrt(kB * temperature / m);

    //--------------------------------------------------
    // inject particles
    //--------------------------------------------------

    for (int p=0; p<num_sim; p++)
    {
        vector3d pos;
        vector3d vel;

        //--------------------------------------------------
        // uniform circular inlet
        //--------------------------------------------------

        double rr =
            inletRadius * sqrt(rnd());

        double theta =
            2.0 * PI * rnd();

        pos.x = rr * cos(theta);
        pos.y = rr * sin(theta);

        pos.z = zInject;

        //--------------------------------------------------
        // thermal velocity
        //--------------------------------------------------

        vel.x =
            sigma * gaussian();

        vel.y =
            sigma * gaussian();

        //--------------------------------------------------
        // IMPORTANT:
        // positive z only
        //--------------------------------------------------

        vel.z =
            sampleHalfMaxwellian(sigma);

        //--------------------------------------------------
        // add particle
        //--------------------------------------------------

        sp.addParticle(pos, vel);
    }
}
/*samples ions at the k=0 edge such that local plasma density equals PLASMA_DEN*/
void Source::injectBeam(int ninject, 
						double beam_radius_mm, 
						double divergence_mrad,  
						double z_inj, 
						double engInitial
						)
{
	/*only sample on global k=0 face*/
	if (domain.mpi_rank!=0) return;
	/** sample particles **/
    double vAbs = sqrt(fabs(2.0*sp.charge*engInitial/sp.mass));
   
	vector3d pos, vel; 

	for (int i=0;i<ninject;i++)
    	{
       		double theta = 2.0 * PI * rnd();

			double r = beam_radius_mm*1e-3 * sqrt(rnd());

			pos.x = r * cos(theta);
			pos.y = r * sin(theta);
			pos.z = z_inj;
		
			
			double sigma = divergence_mrad * 1e-3;

			double div = sigma * gaussian();

			double phi = 2.0*PI*rnd();

			double vr = vAbs * sin(div);

			vel.x = vr*cos(phi);
			vel.y = vr*sin(phi);
			vel.z = vAbs*cos(div);
    
       		sp.addParticle(pos, vel);	
    		
		}    		
}
void Source::sampleCell(
                        int i,
                        int k,
                        double num_real,
                        int spwt,
                        double T
                      )
{
    //--------------------------------------------------
    // mesh spacing
    //--------------------------------------------------

    double dr = domain.dh[0];
    double dz = domain.dh[2];

    //--------------------------------------------------
    // thermal velocity
    //--------------------------------------------------

    double sigma =
        sqrt(kB * T / sp.mass);

    //--------------------------------------------------
    // drift velocity
    //--------------------------------------------------

    vector3d vDrift;

    vDrift.x = 0.0;
    vDrift.y = 0.0;
    vDrift.z = 0.0;

    //--------------------------------------------------
    // number of simulation particles
    //--------------------------------------------------

    int num_sim =
        static_cast<int>(num_real / spwt);

    //--------------------------------------------------
    // stochastic remainder
    //--------------------------------------------------

    double remain =
        num_real/spwt - num_sim;

    if (rnd() < remain)
        num_sim++;

    //--------------------------------------------------
    // cylindrical cell bounds
    //--------------------------------------------------

    double r1 = i * dr;
    double r2 = (i + 1) * dr;

    double z1 = k * dz;
    double z2 = (k + 1) * dz;

    //--------------------------------------------------
    // create particles
    //--------------------------------------------------

    for (int p = 0; p < num_sim; p++)
    {
        vector3d pos;
        vector3d vel;

        //--------------------------------------------------
        // uniform cylindrical sampling
        //--------------------------------------------------

        double rr;

        if (i == 0)
        {
            rr = r2 * sqrt(rnd());
        }
        else
        {
            rr = sqrt(
                        r1*r1 +
                        rnd()*(r2*r2 - r1*r1)
                     );
        }

        //--------------------------------------------------
        // small theta spread
        //--------------------------------------------------
        // helps 2D3V consistency
        //--------------------------------------------------

        double theta =
            (rnd() - 0.5) * 1.0e-3;

        //--------------------------------------------------
        // Cartesian particle position
        //--------------------------------------------------

        pos.x = rr * cos(theta);
        pos.y = rr * sin(theta);

        pos.z = domain.r0[2]+z1 + dz * rnd();

        //--------------------------------------------------
        // Maxwellian velocity
        //--------------------------------------------------

        vel.x =
            sigma * gaussian()
            + vDrift.x;

        vel.y =
            sigma * gaussian()
            + vDrift.y;

        vel.z =
            sigma * gaussian()
            + vDrift.z;

        //--------------------------------------------------
        // add particle
        //--------------------------------------------------

        sp.addParticle(pos, vel);
    }
}



void Source::samplePlasmaChamber_ConstDen(int spwt, double density, double T)
{
    int Nr   = domain.Nr;
    int Nzee = domain.Nz_phys;

    double dr = domain.dh[0];
    double dz = domain.dh[2];
	//--------------------------------------------------
    // drift velocity
    //--------------------------------------------------

    vector3d vDrift;

    vDrift.x = 0.0;
    vDrift.y = 0.0;
    vDrift.z = 0.0;
    //--------------------------------------------------
    // thermal velocity sigma
    //--------------------------------------------------

    double sigma =
		    sqrt(QE * T / sp.mass);

    //--------------------------------------------------
    // loop over RZ cells
    //--------------------------------------------------

    for (int k = 0; k < Nzee-1; k++)
    {
        for (int i = 0; i < Nr; i++)
        {
            //--------------------------------------------------
            // only plasma region
            //--------------------------------------------------

            if (domain.object[i + k*Nr] != OBJ_PLASMA)
                continue;

            //--------------------------------------------------
            // cylindrical cell volume
            //--------------------------------------------------

            double r1;
			double r2;

			if(i==0)
				{
			    r1 = 0.0;
			    r2 = 0.5*dr;
				}
			else
				{
			    r1 = (i-0.5)*dr;
			    r2 = (i+0.5)*dr;
				}
				
			double cellVolume =
			    PI*(r2*r2-r1*r1)*dz;
            //--------------------------------------------------
            // number of REAL particles
            //--------------------------------------------------

            double num_real =
                density * cellVolume;

            //--------------------------------------------------
            // number of simulation particles
            //--------------------------------------------------

            int num_sim =
                static_cast<int>(num_real / spwt);

            //--------------------------------------------------
            // stochastic remainder
            //--------------------------------------------------

            double remain =
                num_real/spwt - num_sim;

            if (rnd() < remain)
                num_sim++;

            //--------------------------------------------------
            // inject particles
            //--------------------------------------------------

            for (int p = 0; p < num_sim; p++)
            {
                //--------------------------------------------------
                // uniform position in annular RZ cell
                //--------------------------------------------------
               double rr =
    				sqrt(
				        r1*r1
        				+
				        rnd()*(r2*r2-r1*r1)
    					);
                double theta =
                    2.0 * PI * rnd();

                vector3d pos;

                pos.x = rr * cos(theta);
                pos.y = rr * sin(theta);
                
                pos.z =
                    	domain.r0[2]+ (k + rnd()) * dz  ;
                    	/** for testing **/ 
		//		if(pos.z < (inputPar.ChamberStart+inputPar.WallThick) 
		//				|| pos.z>(inputPar.ChamberEnd)  )
		//			continue;
				/** for testing **/ 
				//if(pos.z < 25.0e-03 || pos.z>35.0e-03 || rr >8.0e-03 )
				//	continue;	
                //--------------------------------------------------
                // Maxwellian velocity
                //--------------------------------------------------

                vector3d vel;

                vel.x =
                    sigma * gaussian();

                vel.y =
                    sigma * gaussian();

                vel.z =
                    sigma * gaussian();
              vel.z=0.0;
                //--------------------------------------------------
                // optional drift velocity
                //--------------------------------------------------

                vel.x += vDrift.x;
                vel.y += vDrift.y;
                vel.z += vDrift.z;

                //--------------------------------------------------
                // add particle
                //--------------------------------------------------

                sp.addParticle(pos, vel);
            }
        }
    }
}

void initiateIonsPlasmaChamber_ConstDen(
    FluidSpecies& fs,	Domain& domain,
    double density,	double T
	)
{
    int Nr = domain.Nr;
    int Nz = domain.Nz_phys;

    for (int k=0; k<Nz; k++)
    {
        int kg = k + 1;

        for (int i=0; i<Nr; i++)
        {
            if (
                domain.object[i + k*Nr]
                != OBJ_PLASMA
            )
                continue;

            //----------------------------------
            // density
            //----------------------------------

            fs.density(i,kg) =
                density;

            //----------------------------------
            // temperature
            //----------------------------------

            fs.temperature(i,kg) =
                T;

            //----------------------------------
            // pressure
            //----------------------------------

            fs.pressure(i,kg) =
                density * QE * T;

            //----------------------------------
            // zero bulk velocity
            //----------------------------------

            fs.ur(i,kg) = 0.0;

            fs.uz(i,kg) = 0.0;

            //----------------------------------
            // internal energy
            //----------------------------------

            fs.internalEnergy(i,kg) =
                1.5 * density * QE * T;
        }
    }

    //------------------------------------------
    // MPI sync
    //------------------------------------------

    fs.density.updateBoundaries(
        domain.mpi_rank,
        domain.mpi_size,
        1
    );
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



void Source::testScatter(int spwt, double density, double T)
{
    int Nr   = domain.Nr;
    int Nzee = domain.Nz_phys;

    double dr = domain.dh[0];
    double dz = domain.dh[2];
	//--------------------------------------------------
    // drift velocity
    //--------------------------------------------------

    vector3d vDrift;

    vDrift.x = 0.0;
    vDrift.y = 0.0;
    vDrift.z = 0.0;
    //--------------------------------------------------
    // thermal velocity sigma
    //--------------------------------------------------

    double sigma =
		    sqrt(QE * T / sp.mass);

    //--------------------------------------------------
    // loop over RZ cells
    //--------------------------------------------------

    for (int k = 0; k < Nzee-1; k++)
    {
        for (int i = 0; i < Nr-1; i++)
        {
            //--------------------------------------------------
            // only plasma region
            //--------------------------------------------------

            if (domain.object[i + k*Nr] != OBJ_PLASMA)
                continue;

            //--------------------------------------------------
            // cylindrical cell volume
            //--------------------------------------------------

            double r1;
			double r2;

			   
			    r1 = (i+0.65)*dr;
				
            //--------------------------------------------------
            // number of REAL particles

            //--------------------------------------------------
            // inject particles
            //--------------------------------------------------

                double theta =0.0;

                vector3d pos;

                pos.x = r1 * cos(theta);
                pos.y = r1 * sin(theta);
                
                pos.z = domain.r0[2]+ (k + 0.25) * dz  ;
                    	/** for testing **/ 
				if(pos.z < (inputPar.Chamber.ChamberStart
							+inputPar.Chamber.WallThick) 
						|| pos.z>(inputPar.Chamber.ChamberEnd)  )
					continue;
				/** for testing **/ 
				//if(pos.z < 25.0e-03 || pos.z>35.0e-03 || rr >8.0e-03 )
				//	continue;	
                //--------------------------------------------------
                // Maxwellian velocity
                //--------------------------------------------------

                vector3d vel;

                vel.x =0.0;

                vel.y =0.0;
	            vel.z=0.0;
                
                sp.addParticle(pos, vel);
            
        }
    }
}

