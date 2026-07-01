/*defines the simulation domain*/
#ifndef _DOMAIN_H
#define _DOMAIN_H

#include <math.h>
#include <iostream>
#include <sstream>
#include "Field.h"
#include "Utils.h"
//#include "Toroid.h"
/*define constants*/
enum ObjectType
{
    OBJ_VOID        = 0,

    // plasma regions
    OBJ_PLASMA      = 1,
    OBJ_EXTRACTION  = 2,

    // solid materials
    OBJ_WALL        = 11,
    OBJ_ELECTRODE   = 12,
    OBJ_DIELECTRIC  = 31,

    // optional future
    OBJ_INLET       = 40,
    OBJ_FARADAY     = 50
};

class Domain 
{
public:	
	/*VARIABlES*/
	double r0[3];
	double dh[3];
	double rd[3];	/*diagonal corner*/
	int Nr,Ntheta,Nzee;
	int nn;
	int nzee; 		/** nz local **/
	int Nz_global, Nz_phys, Nz_local;
	
	

	Field phi;
	Field rho;
	
	FieldCyl Estatic;      // replaces ef
	FieldCyl Bstatic;      // replaces B

	FieldCyl Erf0;          // RF electric field
	FieldCyl Brf0;          // RF magnetic field

	Field3d Arf;


//	void calcBField_Solenoid(Solenoid &sol);
	
//	double dt;
	double time = 0;
	/*FUNCTIONS*/

	/*constructor, allocats memory*/
	Domain( int Nr, int Ntheta, int Nzee, int Nzlocal);

	/*destructor, frees memory*/
	~Domain();

 	std::vector<ObjectType> object;
	/*MPI*/
	int mpi_rank;
	int mpi_size;
	int mpi_zee;
	int mpi_size_zee;
	int num_nodes_global;
	
	void setOrigin(double r0,double theta0, double zee0) {this->r0[0] = r0;this->r0[1]=theta0;this->r0[2]=zee0;}
	void setSpacing(double dr,double dtheta, double dzee) {this->dh[0]=dr;this->dh[1]=dtheta;this->dh[2]=dzee;}
	/*sets simulation time step*/
//	void setDt(double dt) {this->dt=dt;}
	
	/*converts physical position to logical coordinate*/
	inline void RtoL(double lc[3], double posCyl[3])
	{
		for (int i=0;i<3;i++)
			lc[i] = (posCyl[i]-r0[i])/(dh[i]);
	}

	/*converts physical position to logical coordinate*/
	
	
	inline void Posv3d2lc(double lc[3], vector3d pos)
	{	
	 	vectorCyl pos_c;
		pos_c= Cart2Cyl( pos);
		lc[0] = (pos_c.r-r0[0])/(dh[0]);
		lc[1] = (pos_c.theta-r0[1])/(dh[1]);
		lc[2] = (pos_c.zee-r0[2])/(dh[2]);
	}
	
		/*converts logical coordinate to physical position*/
	void pos(double x[3], double i, double j, double k)
	{
		j=0;
		double lc[3] = {i,j,k};
		pos(x, lc);
	}
	
	/*converts logical coordinate to physical position*/
	void pos(double x[3], double lc[3])
	{
		for (int i=0;i<3;i++)
			x[i] = r0[i]+dh[i]*lc[i];
	}
	void  initMPIDomain(int mpi_rank, int mpi_size);	
	

	inline vectorCyl electricField(
    		int i,
		    int k,
    		double rfFactor
			) 
		{
    	return
        	Estatic(i,k)
      		+ Erf0(i,k)*rfFactor;
		}

	inline vectorCyl magneticField(
		    int i,
    		int k,
   			double rfFactor
			) 
		{
    	return
        	Bstatic(i,k)
		      + Brf0(i,k)*rfFactor;
		}
	
	vectorCyl electricField(
    double lc[3],
    double rfFactor
) ;

	vectorCyl magneticField(
    double lc[3],
    double rfFactor
) ;

	
//protected:
	

};

double rnd();	/*prototype for RNG*/


#endif
