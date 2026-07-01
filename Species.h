#ifndef _SPECIES_H
#define _SPECIES_H

#include <stdio.h>
#include <stdlib.h>

#include <vector>
#include <string>

#include <mpi.h>
#include <omp.h>

#include "Domain.h"
#include "Field.h"

using namespace std;
using namespace Const;
//====================================================
// Particle status
//====================================================

enum ParticleStatus
{
    PARTICLE_LOST        = 0,

    PARTICLE_ACTIVE      = 1,

    PARTICLE_COLLECTED   = 3,

	PARTICLE_DETECTED 	 = 4,
	
    PARTICLE_SEND_LEFT   = 100,
    PARTICLE_SEND_RIGHT  = 200
};


//====================================================
// Particle structure
//====================================================

struct Particle
{
    vector3d pos;
    vector3d vel;

    int STATUS;
};


//====================================================
// Species type
//====================================================

enum SpeciesType
{
    ELECTRON,
    ION,
    NEUTRAL
};


//====================================================
// Base species class
//====================================================

class Species
{
public:

    Species(
        SpeciesType speciesType,
        string name,
        double mass,
        double charge
    )
    :
        speciesType(speciesType),
        name(name),
        mass(mass),
        charge(charge)
    {}

    virtual ~Species() {}

    //------------------------------------------------
    // Common species properties
    //------------------------------------------------

    SpeciesType speciesType;

    string name;

    double mass;
    double charge;
	 double dt;
};



//====================================================
// KINETIC species
//====================================================

class KineticSpecies : public Species
{
public:

    KineticSpecies(
        Domain& domain,
        SpeciesType speciesType,
        string name,
        double mass,
        double charge,
        double spwt
    )
    :
        Species(
            speciesType,
            name,
            mass,
            charge
        ),

        spwt(spwt),

        den(
            domain.Nr,
            domain.nzee
        ),

        domain(domain)

    {
        particles.reserve(100000);
    }

    //------------------------------------------------
    // Particle storage
    //------------------------------------------------

    vector<Particle> particles;

    //------------------------------------------------
    // Macroscopic properties
    //------------------------------------------------

    Field den;

    //------------------------------------------------
    // Particle weight
    //------------------------------------------------

    double spwt;

    //------------------------------------------------
    // Diagnostics
    //------------------------------------------------

    int np      = 0;

    int nLost   = 0;
    int nFC     = 0;

    //------------------------------------------------
    // MPI transfer
    //------------------------------------------------

    int sendL = 0;
    int sendR = 0;

    int recvL = 0;
    int recvR = 0;

    //------------------------------------------------
    // Methods
    //------------------------------------------------

    void addParticle(
        const vector3d& position,
        vector3d velocity
    );

    void checkBounds();

    void applyDomainConditions();

    void transferParticles();

    void computeGasProperties();
    
    void scatterCharge();

protected:

    Domain& domain;

private:

    //------------------------------------------------
    // MPI buffers
    //------------------------------------------------

    Particle* send_bufL = nullptr;
    Particle* send_bufR = nullptr;

    Particle* recv_bufL = nullptr;
    Particle* recv_bufR = nullptr;
};



//====================================================
// FLUID species
//====================================================

class FluidSpecies : public Species
{
public:

    FluidSpecies(
        Domain& domain,
        SpeciesType speciesType,
        string name,
        double mass,
        double charge
    )
    :
        Species(
            speciesType,
            name,
            mass,
            charge
        ),

        density(
            domain.Nr,
            domain.Nz_local
        ),

        pressure(
            domain.Nr,
            domain.Nz_local
        ),

        temperature(
            domain.Nr,
            domain.Nz_local
        ),

        ur(
            domain.Nr,
            domain.Nz_local
        ),

        uz(
            domain.Nr,
            domain.Nz_local
        ),

        internalEnergy(
            domain.Nr,
            domain.Nz_local
        ),
		
		flux_r(
		    domain.Nr,
    		domain.Nz_local
		),

		flux_z(
		    domain.Nr,
    		domain.Nz_local
		),
		density_old(
		    domain.Nr,
    		domain.Nz_local
		),
		source(
		    domain.Nr,
    		domain.Nz_local
		),
        domain(domain)

    {
    	density.clear(); density_old.clear();
		ur.clear();	uz.clear();
		temperature.clear();	pressure.clear();
		internalEnergy.clear();
		flux_r.clear(); 	flux_z.clear();
    }

    //------------------------------------------------
    // Fluid variables
    //------------------------------------------------

    Field density;

    Field pressure;

    Field temperature;

    //------------------------------------------------
    // Velocity
    //------------------------------------------------

    Field ur;
    Field uz;

	//------------------------------------------------
    // flux 
    //------------------------------------------------
	Field flux_r;
    Field flux_z;
    //------------------------------------------------
    // Internal energy
    //------------------------------------------------

    Field internalEnergy;
    
    double dt;


	double mobility=0.05;
	double diffusion;
	
	Field density_old;
	
	Field source;


protected:

	
	
    Domain& domain;
};

#endif