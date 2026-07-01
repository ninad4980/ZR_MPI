#include "MagneticFieldSolver.h"

using namespace Const;

///////////////////////////////////////////////////////////
//
// Constructor
//
///////////////////////////////////////////////////////////

MagneticFieldSolver::
MagneticFieldSolver(
    Domain& domain,
    StaticCoilAssembly& coils
)
:
domain(domain),
coils(coils)
{
}

///////////////////////////////////////////////////////////
//
// Compute magnetic field
//
///////////////////////////////////////////////////////////

void MagneticFieldSolver::compute()
{
    int Nr = domain.Nr;
    int Nz = domain.Nz_phys;

    for(int k=0;k<Nz;k++)
    {
        int kg = k+1;

        for(int i=0;i<Nr;i++)
        {
            //------------------------------------
            // cylindrical mesh location
            //------------------------------------

            vector3d pos;

            pos.x =
                (i+0.5)*
                domain.dh[0];

            pos.y = 0.0;

            pos.z =
                (k+0.5)*
                domain.dh[2];

            //------------------------------------
            // evaluate B
            //------------------------------------

            vector3d B =
                coils.solveB(pos);

            //------------------------------------
            // Cartesian -> cylindrical
            //------------------------------------

            domain.Bstatic(i,kg).r =
                B.x;

            domain.Bstatic(i,kg).theta =
                B.y;

            domain.Bstatic(i,kg).zee =
                B.z;
        }
    }

    domain.Bstatic.updateBoundaries(
        domain.mpi_rank,
        domain.mpi_size,
        1
    );
}