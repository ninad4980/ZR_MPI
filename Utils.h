#ifndef _UTILS_H
#define _UTILS_H


#include <iostream>
#include <iomanip>
#include <random>
#include <fstream>
#include <mpi.h>
#include <type_traits>
#include <math.h>
const double cprecision = 1e-12;

namespace Const
{
	const double C = 3.0e8;					// speed of light
	const double EPS_0 = 8.85418782e-12;  	// C/(V*m), vacuum permittivity
	const double PI = 3.141592653;			// pi
	const double MU_0 = 4.0*PI*1e-7  ;  	// N/A^2 magnetic permeability 
	const double QE = 1.602176565e-19;		// C, electron charge
	const double AMU = 1.660538921e-27;		// kg, atomic mass unit
	const double ME = 9.10938215e-31;		// kg, electron mass
	const double kB = 1.380648e-23;			// J/K, Boltzmann constant
	const double EvToK = QE/kB;				// 1eV in K ~ 11604
	const double EV_TO_K = 11604.52 ;
}
double rnd();

class vector3d
{
public:

    double x,y,z;

    //-----------------------------------
    // constructors
    //-----------------------------------

    vector3d()
    : x(0.0), y(0.0), z(0.0)
    {}

    vector3d(double X,double Y,double Z)
    : x(X), y(Y), z(Z)
    {}

    //-----------------------------------
    // norms
    //-----------------------------------

    double norm() const
    {
        return sqrt(x*x+y*y+z*z);
    }

    double norm2() const
    {
        return x*x+y*y+z*z;
    }

    //-----------------------------------
    // arithmetic
    //-----------------------------------

    vector3d operator+(const vector3d& b) const
    {
        return {x+b.x, y+b.y, z+b.z};
    }

    vector3d operator-(const vector3d& b) const
    {
        return {x-b.x, y-b.y, z-b.z};
    }

    vector3d operator*(double s) const
    {
        return {x*s, y*s, z*s};
    }

    //-----------------------------------
    // dot product
    //-----------------------------------

    double operator*(const vector3d& b) const
    {
        return x*b.x + y*b.y + z*b.z;
    }

    //-----------------------------------
    // cross product
    //-----------------------------------

    vector3d operator||(const vector3d& b) const
    {
        return
        {
            y*b.z - z*b.y,
            z*b.x - x*b.z,
            x*b.y - y*b.x
        };
    }

    //-----------------------------------
    // assignment from scalar
    //-----------------------------------

    vector3d& operator=(double s)
    {
        x=s;
        y=s;
        z=s;
        return *this;
    }
    
    vector3d& operator+=(const vector3d& b)
		{
    	x += b.x;
    	y += b.y;
    	z += b.z;
    	return *this;
		}
	vector3d& operator-=(const vector3d& b)
		{
    	x -= b.x;
    	y -= b.y;
    	z -= b.z;
    	return *this;
		}
	vector3d& operator*=(const double b)
		{
    	x = x * b;
    	y = y * b;
    	z = z * b;
    	return *this;
		}	
};


class vectorCyl
{
public:

    //-----------------------------------
    // cylindrical components
    //-----------------------------------

    double r;
    double theta;
    double zee;

    //-----------------------------------
    // constructors
    //-----------------------------------

    vectorCyl()
    : r(0.0), theta(0.0), zee(0.0)
    {}

    vectorCyl(double R,double T,double Z)
    : r(R), theta(T), zee(Z)
    {}

    //-----------------------------------
    // arithmetic
    //-----------------------------------

    vectorCyl operator+(const vectorCyl& b) const
    {
        return
        {
            r     + b.r,
            theta + b.theta,
            zee   + b.zee
        };
    }

    vectorCyl operator-(const vectorCyl& b) const
    {
        return
        {
            r     - b.r,
            theta - b.theta,
            zee   - b.zee
        };
    }

    vectorCyl operator*(double s) const
    {
        return
        {
            r*s,
            theta*s,
            zee*s
        };
    }

    //-----------------------------------
    // assignment from scalar
    //-----------------------------------

    vectorCyl& operator=(double s)
    {
        r     = s;
        theta = s;
        zee   = s;

        return *this;
    }

    //-----------------------------------
    // inplace operators
    //-----------------------------------

    vectorCyl& operator+=(const vectorCyl& b)
    {
        r     += b.r;
        theta += b.theta;
        zee   += b.zee;

        return *this;
    }

    vectorCyl& operator-=(const vectorCyl& b)
    {
        r     -= b.r;
        theta -= b.theta;
        zee   -= b.zee;

        return *this;
    }
};

struct phasespace
{
    vector3d pos;
    vector3d vel;

    phasespace() = default;

    phasespace(
        const vector3d& p,
        const vector3d& v)
        : pos(p), vel(v)
    {}
};

vector3d Cyl2Cart(vectorCyl X);
vectorCyl Cart2Cyl( vector3d X);
vector3d cylV2cart(vectorCyl A, vectorCyl r);
vectorCyl cartV2cyl(vector3d A,vector3d X);


double gaussian();

double sampleHalfMaxwellian(double sigma);

vector3d sampleMaxwellianVelocity(
    double T,	double mass
);

/** ! Box Muller distribution **/
double BoxMuller();


phasespace BorisPush(
    double dt,
    const vector3d& r,	const vector3d& v,
    const vector3d& E,	const vector3d& B,
    double q,	double m
);

double round_to(double value, double precision );
bool lt (double a, double b);
bool gt (double a, double b);
bool lte (double a, double b);
bool gte (double a, double b);
bool et (double a, double b);
bool net (double a, double b);

#endif




