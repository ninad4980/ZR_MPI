/*Field is a container for mesh node data
It contains functions for scatter/gather, and division by volume*/

#ifndef _FIELD_H
#define _FIELD_H


#include <fstream>
#include <mpi.h>
#include <type_traits>
#include <math.h>
#include <cassert>
#include "Utils.h"
using namespace Const;

template <typename TT>
void foo (TT *arr, int size)
	{
	for (int i=0;i<size;i++)
		{
		arr[i] =0.0;
		}
	}


template <typename T>
class Field_
{
public:
	
	/*constructor*/
	Field_(int nr, int nzee) :
	nr{nr}, nzee{nzee}		/*sets this->domain to domain*/
	{
		data = new T [nr*nzee];
	}
	
	/*destructor, frees memory*/
	~Field_()
	{
		delete[] data;
	}	
	
	
	/* scatters scalar value onto a field at logical coordinate lc*/
	/* scatter scalar value onto RZ field */
	void scatter(double lc[3], T value, double r0, double dr)
	{
    // ----------------------------------------
    // logical coordinates
    // ----------------------------------------
	lc[1]=0.0;

    int i = (int)lc[0];
    double di = lc[0] - i;

    int k = (int)lc[2];
    double dk = lc[2] - k;
    
    int kg = k + 1;

	//compute correction factors
	double ri = r0
				+
				i*dr;
	
	double f1 = (ri+0.5*di*dr)
				/
				(ri+0.5*dr)
				;   
	
	double f2 = (ri+0.5*(di+1)*dr)
				/
				(ri+0.5*dr)
				;
	
	
		data[i + kg*nr] += value*(1-di)*(1.0-dk)*f2;
		
		data[i+1 + kg*nr] += value*(di)*(1-dk)*f1;
		
		data[i + (kg+1)*nr] += value*(1-di)*(dk)*f2;
		
		data[i+1 + (kg+1)*nr] += value*(di)*(dk)*f1;
		
	return;	
			
	}
			


	/* gather field value from RZ field */
	T gather(double lc[3], double r0, double dr)
	{
    lc[1] = 0.0;

    int i = (int)lc[0];
    double di = lc[0] - i;

    int k = (int)lc[2];
    double dk = lc[2] - k;

    int kg = k + 1;
	//--------------------------------------------------	
	//compute correction factors
	//--------------------------------------------------	
		double ri = r0
					+i*dr
					;
		
		double f1 = 
					(ri+0.5*di*dr)
					/
					(ri+0.5*dr)
					;   
		double f2 = (ri+0.5*(di+1)*dr)
					/
					(ri+0.5*dr)
					;

		/*gather electric field onto particle position*/
		
		T val =
      			data[i     + kg*nr]       * (1.0-di)*(1.0-dk)*f2
			    + data[i + 1 + kg*nr]       * di      *(1.0-dk)*f1
    			+ data[i     + (kg+1)*nr]   * (1.0-di)*dk      *f2
    			+ data[i + 1 + (kg+1)*nr]   * di      *dk      *f1
    			;
        		
		    return val;
	}	
	/*shortcut to data*/
	T& at(int i, int k)
	{
    return data[i + k*nr];
	}

	const T& at(int i, int k) const
	{
    return data[i + k*nr];
	}
	
	T& operator()(int i, int k)
	{
    return data[i + k*nr];
		}

	const T& operator()(int i, int k) const
	{
    return data[i + k*nr];
		}

	/*shortcut to data*/
	void clear()
	{
    for (int i=0;i<nr*nzee;i++)
        data[i] = T{};
	}
	
	/*divides value by node volume*/
	void divideByNodeVolume(
    int rank,
    int size,
    double dr,
    double dz)
	{
	
    int Nz_phys = nzee - 2;
		
    	
    for (int k=0;k<Nz_phys;k++)
    	{
    	//--------------------------------------------------
    	// physical nodes only
	    //--------------------------------------------------
    	 int kg = k + 1;
    	 //----------------------------------------------
        // z thickness
        //----------------------------------------------

        double dzNode = dz;

        if (rank > 0 && k == 0)
        	{
            dzNode = 0.5 * dz;
			}
        if (rank < size-1 &&
            	k == Nz_phys-1)
        	{
            dzNode = 0.5 * dz;
        	}
        	
		for (int i=0;i<nr;i++)
			{
				double r1 = std::max(0.0, (i-0.5)*dr);
				double r2 = (i+0.5)*dr;

				double dV =
					    dzNode *
					    Const::PI *
					    (r2*r2-r1*r1);

			if (i==0)
    				dV *= 2.0;		
						
				data[i + kg*nr] /= dV;
			}	
		}			
	}

	/*updates boundaries across MPI processes*/
	void updateBoundaries( int rank, int size, int dim)
	{
	MPI_Status   status;
		
	/** Send to Right **/
	if ( rank < (size-1) ){
        	MPI_Send( &data[(nzee-2)*nr], dim*nr , MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD );
		}
    	/** Rec from left **/
    	if (rank > 0 ){
        	MPI_Recv( &data[0], dim*nr, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &status );
		}	
    	/** Send to left **/
    	if ( rank > 0){
        	MPI_Send( &data[nr], dim*nr, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD );
		}
    	/** Rec from Right **/
    	if ( rank < (size - 1) ){
        	MPI_Recv( &data[(nzee-1)*nr], dim*nr, MPI_DOUBLE,  rank + 1, 1, MPI_COMM_WORLD, &status );
		}
	}
	

	void accumulateBoundaries(int rank,int size,int dim)
	{	
    MPI_Status status;

    if (rank < size-1)
    {
        MPI_Sendrecv(
            &data[(nzee-2)*nr],
            dim*nr,
            MPI_DOUBLE,
            rank+1,
            100,

            &data[(nzee-1)*nr],
            dim*nr,
            MPI_DOUBLE,
            rank+1,
            101,

            MPI_COMM_WORLD,
            &status
        );
    }

    if (rank > 0)
    {
        MPI_Sendrecv(
            &data[nr],
            dim*nr,
            MPI_DOUBLE,
            rank-1,
            101,

            &data[0],
            dim*nr,
            MPI_DOUBLE,
            rank-1,
            100,

            MPI_COMM_WORLD,
            &status
        );
    }

    //--------------------------------------------------
    // add overlapping contributions
    //--------------------------------------------------

    if (rank > 0)
    	{
        for (int i=0;i<dim*nr;i++)
            data[nr+i] += data[i];
    	}

    if (rank < size-1)
    	{
        for (int i=0;i<dim*nr;i++)
            data[(nzee-2)*nr+i] += data[(nzee-1)*nr+i];
    	}
	}
		
	/*performs element by element division by another field*/
	void operator /= (const Field_ &F) {
		for (int i =0; i< nr*nzee; i++)
			{
				if (F.data[i]!=0)
				  data[i] /= F.data[i];
				else
				  data[i] = 0;
			}  
	}
	
	/*performs element by element division by another field*/
	void operator *= (const Field_ &F) {
		for (int i =0; i<nr*nzee; i++)
			{
				
				  data[i] *= F.data[i];

			}  
	}
	
	void multiplyFieldByVal(double val)
	{
	for (int i =0; i<nr*nzee; i++)
			{
				{data[i] = data[i] * val;}
			}
	}
	
	void operator/=(const T d)
	{
	 assert(std::abs(d) > std::numeric_limits<T>::epsilon());
    for (int i = 0; i < nr*nzee; i++)
        data[i] /= d;
	}

	void operator*=(const T d)
	{
    for (int i = 0; i < nr*nzee; i++)
        data[i] *= d;
	}
	void operator+=(const T d)
	{
    for (int i=0; i<nr*nzee; i++)
        data[i] += d;
	}
	Field_ operator/(const T d) const
	{
    Field_ out(*this);

    out /= d;

    return out;
	}

	Field_ operator*(const T d) const
	{
    Field_ out(*this);

    out *= d;

    return out;
	}
	void operator-=(const T d)
	{
    for (int i=0; i<nr*nzee; i++)
        data[i] -= d;
	}
	int nr, ntheta, nzee;

	T *data;	/*data held by this field*/

};

//some typedefs
typedef Field_<double> Field;
typedef Field_<vector3d> Field3d;
typedef Field_<vectorCyl> FieldCyl;

#endif
