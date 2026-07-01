#include <iostream>
#include "Utils.h"

using namespace Const;

/*random number generator*/
std::random_device rd;
std::mt19937 mt_gen(rd());		/*seed*/
std::uniform_real_distribution<double> rnd_dist(0, 1.0);
double rnd()
{
	return rnd_dist(mt_gen);
}

vector3d Cyl2Cart(vectorCyl X)
{

vector3d Y;
double r= X.r; double theta = X.theta; double zee= X.zee;

    Y.x= r*cos(theta) ;
    Y.y= r*sin(theta) ;
    Y.z= zee;

return (Y);
}


vectorCyl Cart2Cyl( vector3d X)
{
vectorCyl Y;
double x= X.x;double y= X.y;double z= X.z;

		Y.r=sqrt(x*x + y*y);

		Y.theta=atan2(y,x);
		if (Y.theta < 0.0){Y.theta = 2.0*PI+Y.theta;}

		Y.zee= X.z;


return(Y);
}


vector3d cylV2cart(vectorCyl A, vectorCyl r)
{
double theta = r.theta; double zee= r.zee;
vector3d B;
	
    B.x = A.r * cos(theta) - A.theta* sin(theta) ;
    B.y = A.r * sin(theta) + A.theta* cos(r.theta) ;
    B.z = A.zee ;

return B;

}

vectorCyl cartV2cyl(vector3d A, vector3d X)
{
vectorCyl B;

	double theta=atan2(X.y , X.x);
		if (theta < 0.0){theta = 2.0*PI+theta;}

    B.r =  A.x * cos(theta)    + A.y* sin(theta) ;
    B.theta = -A.x* sin(theta) + A.y* cos(theta) ;
    B.zee =  A.z ;

return B;

}


phasespace BorisPush(
    double dt,
    const vector3d& r,
    const vector3d& v,
    const vector3d& E,
    const vector3d& B,
    double q,
    double m
)
{
	phasespace ps;

    const double h = 0.5 * dt * q / m;

    //---------------------------------------
    // half electric acceleration
    //---------------------------------------

    vector3d vminus = v + E * h;

    //---------------------------------------
    // magnetic rotation
    //---------------------------------------
    vector3d t = B * h;

    double t2 =
          t.x*t.x
        + t.y*t.y
        + t.z*t.z;

    vector3d s = t * (2.0 / (1.0 + t2));

    vector3d vprime =
        vminus + (vminus || t);

    vector3d vplus =
        vminus + (vprime || s);

    //---------------------------------------
    // second half electric acceleration
    //---------------------------------------

    vector3d vnew =
        vplus + E * h;

    //---------------------------------------
    // advance position
    //---------------------------------------

    vector3d rnew =
        r + vnew * dt;

    ps.pos = rnew;
    ps.vel = vnew;

    return ps;
}

double gaussian()
{
    static bool hasSpare = false;
    static double spare;

    if(hasSpare)
    {
        hasSpare = false;
        return spare;
    }

    hasSpare = true;

    double u,v,s;

    do
    {
        u = 2.0*rnd() - 1.0;
        v = 2.0*rnd() - 1.0;
        s = u*u + v*v;

    } while(s >= 1.0 || s == 0.0);

    s = sqrt(-2.0*log(s)/s);

    spare = v*s;

    return u*s;
}


vector3d sampleMaxwellianVelocity(
    double T,
    double mass
)
{
    double sigma =
        sqrt(
            Const::kB*T/mass
        );

    vector3d v;

    v.x = sigma * BoxMuller();
    v.y = sigma * BoxMuller();
    v.z = sigma * BoxMuller();

    return v;
}


double sampleHalfMaxwellian(double sigma)
{
    double v;

    do
    {
        v = sigma * gaussian();
    }
    while (v <= 0.0);

    return v;
}


/** ! Box Muller distribution **/
double BoxMuller()
{
    double r1 = rnd();
    double r2 = rnd();

    return sqrt(-2.0*log(r1))
           * cos(2.0*Const::PI*r2);
}


double round_to(double value, double precision )
{
    return std::round(value / precision) * precision;
}

bool lt (double a, double b)
{
    if( round_to(a, cprecision) <  round_to(b, cprecision)  )
    return true;
    else return false;
}

bool gt (double a, double b)
{
    if( round_to(a, cprecision) >  round_to(b, cprecision)  )
    return true;
    else return false;
}

bool lte (double a, double b)
{
    if( round_to(a, cprecision) <=  round_to(b, cprecision)  )
    return true;
    else return false;
}

bool gte (double a, double b)
{
    if( round_to(a, cprecision) >=  round_to(b, cprecision)  )
    return true;
    else return false;
}

bool et (double a, double b)
{
    if( round_to(a, cprecision) ==  round_to(b, cprecision)  )
    return true;
    else return false;
}

bool net (double a, double b)
{
    if( round_to(a, cprecision) !=  round_to(b, cprecision)  )
    return true;
    else return false;
}