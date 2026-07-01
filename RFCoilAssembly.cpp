
#include <fstream>
#include "RFCoilAssembly.h"
#include "InputPar.h"
using namespace std;
using namespace Const;


RFCoilAssembly::RFCoilAssembly(const InputPar& input)
{
    radius     = input.RF.Radius;
    turns      = input.RF.Turns;
    pitch	   = input.RF.Pitch;	
    layers     = input.RF.Layers;
    frequency  = input.RF.Frequency;
    power      = input.RF.Power;
    phase      = input.RF.Phase;
    
    
    current =	input.RF.Current;
    
    center =
		{
		    0,
    		0,
		    input.RF.StartZ
		};

    normal = {0,0,1};
    
      build();
      
}

void RFCoilAssembly::build()
{
    loops.clear();
	for(int l=0;l<layers;l++)
    	{
    	double r =
	    	radius+l*pitch;
	    
	    for(int n=0;n<turns;n++)
		    {
    	    CurrentLoop loop;

    	    loop.radius = radius;

    	    loop.current = current;

        	loop.normal = normal;

	        loop.center = center
     	               +
        	              normal
            	        *
                	      (
                    	   (n-(turns-1)/2.0)
                       			*
		                       	pitch
                      		);

        		loops.push_back(loop);
    			}
    	}		
}

vector3d RFCoilAssembly::solveB(
const vector3d& pos
) const
{
    vector3d B;

    B.x=0;
    B.y=0;
    B.z=0;

    for(const auto& loop:loops)
        B += loop.solveB(pos);

    return B;
}

vector3d RFCoilAssembly::solveA(
const vector3d& pos
) const
{
    vector3d A;

    A.x=0;
    A.y=0;
    A.z=0;

    for(const auto& loop:loops)
        A += loop.solveA(pos);

    return A;
}

