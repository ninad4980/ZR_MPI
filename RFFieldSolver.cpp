#include "RFFieldSolver.h"

using namespace Const;

RFFieldSolver::RFFieldSolver(
    Domain& domain,
    RFCoilAssembly& antenna
)
:
domain(domain),
antenna(antenna)
{
}

///////////////////////////////////////////////////////////
//
// Compute RF vector potential
//
///////////////////////////////////////////////////////////

void RFFieldSolver::computeA()
{
    int Nr = domain.Nr;
    int Nz = domain.Nz_phys;

    for(int k=0;k<Nz;k++)
    {
        int kg = k+1;

        for(int i=0;i<Nr;i++)
        {
            //-------------------------------------------------
            // cell center
            //-------------------------------------------------

            vector3d pos;

            pos.x = (i+0.5)*domain.dh[0];
            pos.y = 0.0;
            pos.z = (k+0.5)*domain.dh[2];

            //-------------------------------------------------
            // RF vector potential
            //-------------------------------------------------

            vector3d A =
                antenna.solveA(pos);

            domain.Arf(i,kg).x = A.x;
            domain.Arf(i,kg).y = A.y;
            domain.Arf(i,kg).z = A.z;
        }
    }

    domain.Arf.updateBoundaries(
        domain.mpi_rank,
        domain.mpi_size,
        1
    );
}

///////////////////////////////////////////////////////////
//
// Compute RF magnetic field
//
///////////////////////////////////////////////////////////

void RFFieldSolver::computeB()
{
    int Nr = domain.Nr;
    int Nz = domain.Nz_phys;

    for(int k=0;k<Nz;k++)
    {
        int kg = k+1;

        for(int i=0;i<Nr;i++)
        {
            vector3d pos;

            pos.x=(i+0.5)*domain.dh[0];
            pos.y=0.0;
            pos.z=(k+0.5)*domain.dh[2];

			vector3d B =
			    antenna.solveB(pos);	
			    
         	/********  y=0 assumed *****/
           	domain.Brf0(i,kg).r      = B.x;
			domain.Brf0(i,kg).theta  = B.y;
			domain.Brf0(i,kg).zee    = B.z;
        }
    }

    domain.Brf0.updateBoundaries(
        domain.mpi_rank,
        domain.mpi_size,
        1
    );
}

///////////////////////////////////////////////////////////
//
// RF electric field
//
///////////////////////////////////////////////////////////

void RFFieldSolver::computeERF()
{
    const double omega =
        2.0*
        PI*
        antenna.frequency;

    int Nr = domain.Nr;
    int Nz = domain.Nz_phys;

    for(int k=0;k<Nz;k++)
    {
        int kg=k+1;

        for(int i=0;i<Nr;i++)
        {
            domain.Erf0(i,kg).r=0.0;

    		domain.Erf0(i,kg).theta=
        				omega*
				        domain.Arf(i,kg).y;

		    domain.Erf0(i,kg).zee=0.0;
        }
    }

    domain.Erf0.updateBoundaries(
        domain.mpi_rank,
        domain.mpi_size,
        1
    );
}



/*
RFFieldSolver::update(double time)
{
	const double omega =
        2.0*
        Const::PI*
        antenna.frequency;

    const double factor =
        std::cos(
            omega*time
        );

    int Nr = domain.Nr;
    int Nz = domain.Nz_phys;

    for(int k=0;k<Nz;k++)
    {
        int kg = k+1;

        for(int i=0;i<Nr;i++)
        {
            domain.Erf(i,kg).r =
                factor*
                domain.Erf(i,kg).r;

            domain.Erf(i,kg).theta =
                factor*
                domain.Erf(i,kg).theta;

            domain.Erf(i,kg).zee =
                factor*
                domain.Erf(i,kg).zee;

            domain.Brf(i,kg).r =
                factor*
                domain.Brf(i,kg).r;

            domain.Brf(i,kg).theta =
                factor*
                domain.Brf(i,kg).theta;

            domain.Brf(i,kg).zee =
                factor*
                domain.Brf(i,kg).zee;
        }
    }

}
*/