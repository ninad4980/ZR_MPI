#include "PotentialSolver.h"

#include <math.h>
#include <iostream>

/**
	0 = outside domain
	1 = inside domain
	
	10 = central nodes
	
	11 = r max
	12 = K min / kmax

	21 = fixed potential
			 
**/	
using namespace std;
void PotentialSolver::calcRhoC(
    std::vector<std::unique_ptr<Species>> &speciesList
)
{
    Field &rho = domain.rho;

    rho.clear();

    int N =
        domain.Nr *
        domain.Nz_local;

    //--------------------------------------------------
    // loop over all species
    //--------------------------------------------------

    for (auto &speciesPtr : speciesList)
    {
        Species *sp = speciesPtr.get();

        //--------------------------------------------------
        // skip neutrals
        //--------------------------------------------------

        if (sp->charge == 0.0)
            continue;

        //--------------------------------------------------
        // kinetic species
        //--------------------------------------------------

        if (
            KineticSpecies *ks =
            dynamic_cast<KineticSpecies*>(sp)
        )
        {
            for (int i=0; i<N; i++)
            {
                rho.data[i] +=
                    ks->charge *
                    ks->den.data[i];
            }
        }

        //--------------------------------------------------
        // fluid species
        //--------------------------------------------------

        else if (
            FluidSpecies *fs =
            dynamic_cast<FluidSpecies*>(sp)
        )
        {
            for (int i=0; i<N; i++)
            {
                rho.data[i] +=
                    fs->charge *
                    fs->density.data[i];
            }
        }
    }

    //--------------------------------------------------
    // synchronize MPI ghost nodes
    //--------------------------------------------------

    rho.updateBoundaries(
        domain.mpi_rank,
        domain.mpi_size,
        1
    );
}

/*solves Poisson equation using Gauss-Seidel*/
bool PotentialSolver::solvePhi()
{
    Field &phi = domain.phi;
    Field &rho = domain.rho;
	
    const int Nr   = domain.Nr;
    const int Nz   = domain.Nz_local;

    const double dr = domain.dh[0];
    const double dz = domain.dh[2];

    const double invdr2 = 1.0/(dr*dr);
    const double invdz2 = 1.0/(dz*dz);

    const double eps0 = EPS_0;

    const int rank = domain.mpi_rank;
    const int size = domain.mpi_size;

    w               = 1.4;
    tolerance       = 1.0e-8;
    max_solver_it   = 10000;

    converged = false;

    // temporary old solution
    std::vector<double> phiOld(Nr*Nz,0.0);

    for (int it=0; it<max_solver_it; it++)
    {
        //------------------------------------------------------
        // Exchange ghost cells
        //------------------------------------------------------
        phi.updateBoundaries(rank,size,1);

        //------------------------------------------------------
        // Store old potential
        //------------------------------------------------------
        for (int k=0;k<Nz;k++)
            for (int i=0;i<Nr;i++)
                phiOld[i + k*Nr] = phi.at(i,k);

        //------------------------------------------------------
        // Interior solve
        //------------------------------------------------------
        for (int k=1; k<Nz-1; k++)
        {
            for (int i=0; i<Nr; i++)
            {
                int obj = domain.object[i + (k-1)*Nr];

                //--------------------------------------------------
                // Plasma region
                //--------------------------------------------------
                if (obj == OBJ_PLASMA ||
					    obj == OBJ_EXTRACTION
					    )
                {
                    double phiNew;

                    //--------------------------------------------------
                    // Axis r = 0
                    //--------------------------------------------------
                    if (i == 0)
                    {
                        double A =
                            2.0*invdr2 +
                            2.0*invdz2;

                        phiNew =
                        (
                            2.0*invdr2 * phi.at(1,k)
                            +
                            invdz2 * phi.at(0,k+1)
                            +
                            invdz2 * phi.at(0,k-1)
                            +
                            rho.at(0,k)/eps0
                        ) / A;
                    }

                    //--------------------------------------------------
                    // Outer radial boundary
                    //--------------------------------------------------
                    else if (i == Nr-1)
                    {
                        // Neumann wall
                        phiNew = phi.at(i-1,k);
                    }

                    //--------------------------------------------------
                    // Interior cylindrical nodes
                    //--------------------------------------------------
                    else
                    {
                        double r = i*dr;

                        double Arp =
                            invdr2 + 1.0/(2.0*r*dr);

                        double Arm =
                            invdr2 - 1.0/(2.0*r*dr);

                        double A =
                            2.0*invdr2 +
                            2.0*invdz2;

                        phiNew =
                        (
                            Arp * phi.at(i+1,k)
                            +
                            Arm * phi.at(i-1,k)
                            +
                            invdz2 * phi.at(i,k+1)
                            +
                            invdz2 * phi.at(i,k-1)
                            +
                            rho.at(i,k)/eps0
                        ) / A;
                    }

                    //--------------------------------------------------
                    // SOR update
                    //--------------------------------------------------
                    phi(i,k) =
                        phi.at(i,k)
                        +
                        w*(phiNew - phi.at(i,k));
                }

                //--------------------------------------------------
                // Extraction electrode
                //--------------------------------------------------
                else if (obj == 11)
                {
                    phi(i,k) = phi.at(i,k);
                }

                //--------------------------------------------------
                // Grounded wall
                //--------------------------------------------------
                else if (obj == 10)
                {
                    phi(i,k) = phi.at(i,k);
                }

                //--------------------------------------------------
                // Floating / insulator
                //--------------------------------------------------
                else if (obj == 31)
                {
                    // zero normal gradient
                    if (i > 0)
                        phi(i,k) = phi.at(i-1,k);
                }
            }
        }

        //------------------------------------------------------
        // Axial ghost boundaries
        //------------------------------------------------------

        // left physical boundary
        if (rank == 0)
        {
            int k = 0;

            for (int i=0;i<Nr;i++)
            {
                phi(i,k) = phi.at(i,1);
            }
        }

        // right physical boundary
        if (rank == size-1)
        {
            int k = Nz-1;

            for (int i=0;i<Nr;i++)
            {
                phi(i,k) = phi.at(i,Nz-2);
            }
        }

        //------------------------------------------------------
        // Convergence check
        //------------------------------------------------------
        double localMax = 0.0;

        for (int k=1;k<Nz-1;k++)
        {
            for (int i=0;i<Nr;i++)
            {
                double err =
                    fabs(
                        phi.at(i,k)
                        -
                        phiOld[i + k*Nr]
                    );

                if (err > localMax)
                    localMax = err;
            }
        }

        double globalMax = 0.0;

        MPI_Allreduce(
            &localMax,
            &globalMax,
            1,
            MPI_DOUBLE,
            MPI_MAX,
            MPI_COMM_WORLD
        );
		/*
        if (rank==0 && it%100==0)
        {
            std::cout
                << "Iteration "
                << it
                << " Error = "
                << globalMax
                << "\n";
        }
		*/
        if (globalMax < tolerance)
        {
            converged = true;
			/*
            if (rank==0)
            {
                std::cout
                    << "Poisson converged in "
                    << it
                    << " iterations\n";
            }
			*/
            break;
        }
    }

    //----------------------------------------------------------
    // Final ghost update
    //----------------------------------------------------------
    phi.updateBoundaries(rank,size,1);

    return converged;
}

	
/*computes electric field = -gradient(phi) using 2nd order differencing*/
void PotentialSolver::computeEF()
{
    const int Nr = domain.Nr;

    // physical local planes
    const int Nz = domain.Nz_phys;

    const double dr = domain.dh[0];
    const double dz = domain.dh[2];

    Field& phi = domain.phi;

    FieldCyl& ef = domain.Estatic;

    ef.clear();

    //--------------------------------------------------
    // Compute electric field
    //--------------------------------------------------

    for (int k=0; k<Nz; k++)
    {
        // field index including ghost offset
        int kg = k + 1;

        for (int i=0; i<Nr; i++)
        {
        	int obj =
                domain.object[i + k*Nr];
            //--------------------------------------------------
            // plasma region only
            //--------------------------------------------------

            if (obj != 1 && obj != 2)
                continue;

            //--------------------------------------------------
            // Er
            //--------------------------------------------------

            if (i == 0)
            {
                //--------------------------------------------------
                // axis symmetry:
                // dphi/dr = 0 at r=0
                //--------------------------------------------------

                ef(0,kg).r = 0.0;
            }
            else if (i == Nr-1)
            {
                //--------------------------------------------------
                // outer radial boundary
                // backward difference
                //--------------------------------------------------

                ef(i,kg).r =
                    -(phi.at(i,kg)
                    - phi.at(i-1,kg))
                    /dr;
                    
            }
            else
            {
                //--------------------------------------------------
                // centered difference
                //--------------------------------------------------

                ef(i,kg).r =
                    -(phi.at(i+1,kg)
                    - phi.at(i-1,kg))
                    /(2.0*dr); 
            }

            //--------------------------------------------------
            // Etheta = 0 in axisymmetric RZ
            //--------------------------------------------------

            ef(i,kg).theta = 0.0;

            //--------------------------------------------------
            // Ez
            //--------------------------------------------------

            //--------------------------------------------------
            // left physical boundary
            //--------------------------------------------------

            if (k == 0)
            {
                ef(i,kg).zee =
                    -(phi.at(i,kg+1)
                    - phi.at(i,kg))
                    /dz;

            }

            //--------------------------------------------------
            // right physical boundary
            //--------------------------------------------------

            else if (k == Nz-1)
            {
                ef(i,kg).zee =
                    -(phi.at(i,kg)
                    - phi.at(i,kg-1))
                    /dz;
            }

            //--------------------------------------------------
            // interior
            //--------------------------------------------------

            else
            {
                ef(i,kg).zee =
                    -(phi.at(i,kg+1)
                    - phi.at(i,kg-1))
                    /(2.0*dz);
            }
        }
    }

    //--------------------------------------------------
    // Exchange ghost cells
    //--------------------------------------------------

    ef.updateBoundaries(
        domain.mpi_rank,
        domain.mpi_size,
        3
    );
}


