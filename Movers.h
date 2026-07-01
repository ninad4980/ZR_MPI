#ifndef _MOVERS_H
#define _MOVERS_H

#include "Species.h"
#include "Domain.h"

#include <memory>

//====================================================
// Base mover class
//====================================================

class Mover
{
public:

    virtual void move(double rfFactor) = 0;

    virtual ~Mover() {}
};


//====================================================
// PIC mover
// electrons
//====================================================

class PICMover : public Mover
{
public:

    PICMover(
        KineticSpecies& species,
        Domain& domain
    )
    :
        sp(species),
        domain(domain)
    {}

    void move(double rfFactor) override;

private:

    KineticSpecies& sp;

    Domain& domain;
};


//====================================================
// DSMC mover
// neutrals
//====================================================

class DSMCMover : public Mover
{
public:

    DSMCMover(
        KineticSpecies& species,
        Domain& domain
    )
    :
        sp(species),
        domain(domain)
    {}

    void move(double) override;

private:

    KineticSpecies& sp;

    Domain& domain;

    //------------------------------------------------
    // DSMC collision routine
    //------------------------------------------------

    void performCollisions();
};


//====================================================
// Fluid mover
// ions
//====================================================

class FluidMover : public Mover
{
public:

    FluidMover(
        FluidSpecies& species,
        Domain& domain
    )
    :
        fs(species),
        domain(domain)
    {}

    void move(double rfFactor) override;
    
    	void computeIonizationSource(
    const KineticSpecies& electrons,
    const KineticSpecies& neutrals
		);

private:

	FluidSpecies& fs;
 
    Domain& domain;

    //------------------------------------------------
    // Fluid update steps
    //------------------------------------------------
    
    void computeVelocity(double rfFactor);
	
	void computeFluxes();

	void advanceDensity();



//    void updateMomentum();

//    void updateEnergy();
};

#endif