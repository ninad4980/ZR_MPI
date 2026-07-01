
#include <fstream>
#include "CurrentLoop.h"
using namespace std;
using namespace Const;

CurrentLoop::CurrentLoop()
{
    center = {0,0,0};

    normal = {0,0,1};

    radius = 0.01;

    current = 1.0;

    integrationSteps = 100;
}


void CurrentLoop::buildLocalBasis(
    vector3d& e1,
    vector3d& e2
) const
{
    if(std::abs(normal.x)<1e-12 &&
       std::abs(normal.y)<1e-12)
    {
        e1 = {0,1,0};
    }
    else
    {
        double s =
        std::sqrt(
            normal.x*normal.x+
            normal.y*normal.y
        );

        e1 =
        {
            -normal.y/s,
             normal.x/s,
             0
        };
    }

    e2 = normal || e1;

	double n = e2.norm();

	if(n>1e-15)
    	e2 *= 1.0/n;
    
}
/** Calculation of magnetic field B-vector for defined position
 and coil: input = position, output = Bfield */
vector3d CurrentLoop::solveB(
const vector3d& position
) const      
{
		vector3d e1,e2;

		buildLocalBasis(e1,e2);
        
        vector3d B,rv,re,dl;        
	
       const double coeff =
					Const::MU_0*
						current/	
						(
						4.0*
						Const::PI
						);
	
        B.x=0.0;
        B.y=0.0;
        B.z=0.0;
        
		double phi;

		const double dphi =
			2.0*Const::PI/
				integrationSteps;
        
        for(int i=0;i<integrationSteps;i++)
        {

                phi=i*dphi;

                re=center+(e1*cos(phi)+e2*sin(phi))*radius;

                rv=position-re;

				double r = rv.norm();

				if(r<1e-12)
				    continue;

				rv *= 1.0/r;
               // dl vector is simple derivation of re-vector in respet to phi 
	       		
				double r3 = r*r*r;
	       dl =
				(
			    e2*cos(phi)
			    -
			    e1*sin(phi)
				)
				*
				(radius*dphi);
                    
                B +=  (dl||rv) * (coeff/r3);                 
          }
	
       return B;
}

vector3d CurrentLoop::solveA(
const vector3d& position
) const      
{	

		vector3d e1,e2;

		buildLocalBasis(e1,e2);

        vector3d A,rv,re,dl;        

        const double coeff =
					Const::MU_0*
						current/	
						(
						4.0*
						Const::PI
						);
        
	A.x=0.0;
	A.y=0.0;
	A.z=0.0;
	
       double phi;

		const double dphi =
				2.0*Const::PI/
					integrationSteps;
        
        for(int i=0;i<integrationSteps;i++)
        {

                phi=i*dphi;

                re=center+(e1*cos(phi)+e2*sin(phi))*radius;

                rv=position-re;
                double r = rv.norm();

				if(r<1e-12)
				    continue;

				rv *= 1.0/r;

               // dl vector is simple derivation of re-vector in respet to phi 
	       
	         dl =
				(
			    e2*cos(phi)
			    -
			    e1*sin(phi)
				)
				*
				(radius*dphi);
                    
	     	A +=  dl * (coeff /r);                   
          }
	
       return A;
}
