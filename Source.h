/*
Definition for a particle source
 */
#ifndef _SOURCE_H
#define _SOURCE_H

#include "Utils.h"
#include "Field.h"
#include "Species.h"
#include "Domain.h"
#include "InputPar.h"
#include <math.h>

using namespace std;


class Source 
{
public:

	
	KineticSpecies &sp;    /*injection species*/

    Source(Domain &domain, KineticSpecies &sp):
		domain(domain), sp(sp)
    	{

    	}
    
	void injectBeam(int ninject, double beam_radius, double beamDiv,  
						double zInitial, double engInitial);
						
	//--------------------------------------------------
	// DSMC inlet injection
	//--------------------------------------------------

	void injectNeutralInlet(
		    double pressure,
    		double temperature,
    		double inletRadius,
    		double zInject,
   			double dt,
   			double massFlow = -1.0
		);					
	
	
	void sampleVolume(int ninject,double radiusMax,
    						double T);
    						
    void samplePlasmaChamber_ConstDen(int spwt,double density, double T);
						
	void sampleCell(int i, int k, double num_real, int spwt, double T );
	
//	void initiateIonsPlasmaChamber_ConstDen( double density0,
//			    double Ti_eV);
	void testScatter(int spwt, double density, double T);

protected:

	Domain &domain;

};

////////////////////////////////////////////////////////
// FREE FUNCTION
////////////////////////////////////////////////////////

void initiateIonsPlasmaChamber_ConstDen(
    FluidSpecies& fs,
    Domain& domain,
    double density,
    double T
);

#endif


