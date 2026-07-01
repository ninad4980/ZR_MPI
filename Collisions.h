#ifndef _COLLISIONS_H
#define _COLLISIONS_H

#include <memory>
#include <vector>
#include "Utils.h"
#include "Species.h"
#include "Domain.h"

//////////////////////////////////////////////////////
// Scatter velocity
//////////////////////////////////////////////////////
// Simple Scatter 
//////////////////////////////////////////////////////
vector3d randomScatter(
		const vector3d& vin
		);
// with thermanl Scatter 
vector3d randomScatter2(
		const vector3d& ve,
	    double Tn,
		double mn
		);
//////////////////////////////////////////////////////


class Interaction
{
public:

    virtual void apply(double dt)=0;

    virtual ~Interaction(){}
};

class CollisionManager
{
public:

    std::vector<
        std::unique_ptr<Interaction>
    > interactions;

    void apply(double dt)
    {
        for(auto& c : interactions)
            c->apply(dt);
    }
};

//////////////////////////////////////////////////////
// Electron-neutral MCC
//////////////////////////////////////////////////////

class ElectronNeutralMCC
    : public Interaction
{
public:

    ElectronNeutralMCC(
        KineticSpecies& electrons,
        KineticSpecies& neutrals,
        FluidSpecies& ions,
        Domain& domain
    );

    void apply(double dt);
	
	vector3d sampleIonizationElectron();
	double sigmaExcitation(double E);

private:

    KineticSpecies& electrons;
    KineticSpecies& neutrals;
    FluidSpecies& ions;

    Domain& domain;

    double sigmaElastic(double EeV);

    double sigmaIonization(double EeV);
};

//////////////////////////////////////////////////////
// Ion-neutral charge exchange
//////////////////////////////////////////////////////

class MCC_CEX
    : public Interaction
{
public:

    MCC_CEX(
        FluidSpecies& ions,
        KineticSpecies& neutrals,
        Domain& domain
    );

    void apply(double dt);

private:

    FluidSpecies& ions;
    KineticSpecies& neutrals;

    Domain& domain;
};

//////////////////////////////////////////////////////
// Neutral-neutral DSMC
//////////////////////////////////////////////////////

class DSMC_MEX : public Interaction
{
public:

    DSMC_MEX(
        KineticSpecies& species,
        Domain& domain
    );


    void apply(double dt);

private:

    void collide(
        vector3d& vel1,
        vector3d& vel2,
        double mass1,
        double mass2
    );

    double evalSigma(
        double g_rel
    );

private:

    double sigma_cr_max = 1e-14;

    double mr;

    double c[4];

    KineticSpecies& species;

    Domain& domain;
};



#endif