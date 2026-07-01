/*defines the simulation domain*/
#include "Domain.h"

#include <math.h>
#include <string>	

using namespace Const;

/*constructor*/
Domain::Domain(int Nr, int Ntheta, int Nzee, int nzee) :
	Nr{Nr}, Ntheta{Ntheta}, Nzee{Nzee}, nzee{nzee},
	phi(Nr,nzee),rho(Nr,nzee),  Estatic(Nr,nzee), Bstatic(Nr,nzee),
	Erf0(Nr,nzee), Brf0(Nr,nzee), Arf(Nr,nzee)
{
	
	phi.clear();
	rho.clear();
	
	Estatic.clear();
	Bstatic.clear();

	Erf0.clear();          // RF electric field
	Brf0.clear();          // RF magnetic field

	Arf.clear();
}
/*destructor, frees memory*/
Domain::~Domain()
{

}

void Domain:: initMPIDomain(int mpi_rank, int mpi_size)
	{
		

		this->mpi_rank = mpi_rank;
	    this->mpi_size = mpi_size;
		/*now set our position in the MPI world*/
		mpi_zee = mpi_rank;

		 //----------------------------------------
    	// save global values
    	//----------------------------------------

	    Nz_global = Nzee;
		/*retain global nodes number*/
		/*this is for clarity*/
		int Nr_global = Nr;
		int Nzee_global = Nzee;

		//----------------------------------------
    	// local physical nodes
    	//----------------------------------------

    	Nz_phys =
        	(Nz_global - 1)/mpi_size + 1;

	    //----------------------------------------
	    // local INCLUDING ghost cells
	    //----------------------------------------

	    Nz_local = Nz_phys + 2;

	    //----------------------------------------
	    // local extents
	    //----------------------------------------
		r0[0] =0.0;
		r0[1] =0.0;	
    	r0[2] = mpi_rank * (Nz_phys - 1) * dh[2];
		rd[0] = r0[0] + (Nr-1) *dh[0];
		rd[1] = r0[1] ;
    	rd[2] = r0[2] + (Nz_phys-1)*dh[2];

		std :: cout<<mpi_rank<<"=rank max \t"<<rd[0]<<"\t"<<rd[1]*180./PI<<"\t"<<rd[2]<<"\n";
		std :: cout<<mpi_rank<<"=rank min \t"<<r0[0]<<"\t"<<r0[1]*180./PI<<"\t"<<r0[2]<<"\n";
		std :: cout<<mpi_rank<<"=rank dh \t"<<dh[0]<<"\t"<<dh[1]*180./PI<<"\t"<<dh[2]<<"\n";
		
		/*object nodes*/
		object.resize(Nr * Nz_phys);

		for (int i=0; i<Nr*Nz_phys; i++)
			{
		    object[i] = OBJ_VOID;
			}
	}

vectorCyl Domain::electricField(
    double lc[3],
    double rfFactor
) 
{
    vectorCyl Es =
        Estatic.gather(
            lc,
            0.0,
            dh[0]
        );

    vectorCyl Erf =
        Erf0.gather(
            lc,
            0.0,
            dh[0]
        );

    return Es + Erf*rfFactor;
}

vectorCyl Domain::magneticField(
    double lc[3],
    double rfFactor
) 
{
    vectorCyl Bs =
        Bstatic.gather(
            lc,
            0.0,
            dh[0]
        );

    vectorCyl Brf =
        Brf0.gather(
            lc,
            0.0,
            dh[0]
        );

    return Bs + Brf*rfFactor;
}

