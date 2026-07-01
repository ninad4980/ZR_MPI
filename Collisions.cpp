#include <vector>
#include <math.h>
#include <iostream>
#include "Collisions.h"
using namespace std;


MCC_CEX :: MCC_CEX(
    FluidSpecies& ions, KineticSpecies& neutrals,	Domain& domain
	) 	:
		ions(ions),	neutrals(neutrals),	domain(domain)
{
}

DSMC_MEX::DSMC_MEX(
    KineticSpecies& species,	Domain& domain
	)	:
		species(species),	domain(domain)
{
    mr =
        species.mass * species.mass
        /
        (2.0 * species.mass);

    c[0] = 4.07e-10;
    c[1] = 0.77;

    c[2] =
        2.0 *
        Const::kB *
        273.15
        /
        mr;

    c[3] =
        std::tgamma(
            2.5 - c[1]
        );
}

ElectronNeutralMCC :: ElectronNeutralMCC(
	    KineticSpecies& electrons,	KineticSpecies& neutrals,	FluidSpecies& ions,
        Domain& domain
		)	:
			electrons(electrons),	neutrals(neutrals),	ions(ions),
			domain(domain)
{
}
 
vector3d ElectronNeutralMCC::
sampleIonizationElectron()
{
    double Esec = 2.0; // eV

    double v =
        sqrt(
            2.0*Const::QE*Esec
            /
            electrons.mass
        );

    double mu =
        2.0*rnd()-1.0;

    double phi =
        2.0*Const::PI*rnd();

    double sint =
        sqrt(1.0-mu*mu);

    vector3d vel;

    vel.x = v*sint*cos(phi);
    vel.y = v*sint*sin(phi);
    vel.z = v*mu;

    return vel;
}

double ElectronNeutralMCC::sigmaExcitation(double E)
{
    if(E < 6.0)
        return 0.0;

    return 1.5e-20*(1.0-6.0/E);
}
double ElectronNeutralMCC	::	sigmaElastic(double E)
{
    return 5e-20;
}



double ElectronNeutralMCC:: sigmaIonization(double E)
{
    if(E < 15.8)
        return 0.0;

    return 2e-20*
           (1.0 - 15.8/E);
}


void ElectronNeutralMCC:: apply(double dt)
{
    double lc[3];

    //--------------------------------------------------
    // maximum collision frequency
    //--------------------------------------------------

    const double nuMax = 5.0e8;

    double Pmax =
        1.0-exp(-nuMax*dt);

    std::vector<Particle> newElectrons;

    for(auto& ep : electrons.particles)
    {
        if(ep.STATUS != PARTICLE_ACTIVE)
            continue;

        //----------------------------------------------
        // null collision selection
        //----------------------------------------------

        if(rnd() > Pmax)
            continue;

        //----------------------------------------------
        // local neutral density
        //----------------------------------------------

        domain.Posv3d2lc(
            lc,
            ep.pos
        );

        double nn =
            neutrals.den.gather(
                lc,
                0.0,
                domain.dh[0]
            );

        //----------------------------------------------
        // electron energy
        //----------------------------------------------

        double v =
            ep.vel.norm();

        double Ee =
            0.5*
            electrons.mass*
            v*v/
            Const::QE;

        //----------------------------------------------
        // cross sections
        //----------------------------------------------

        double sig_el =
            sigmaElastic(Ee);

        double sig_ex =
            sigmaExcitation(Ee);

        double sig_ion =
            sigmaIonization(Ee);

        //----------------------------------------------
        // real frequencies
        //----------------------------------------------

        double nu_el =
            nn*sig_el*v;

        double nu_ex =
            nn*sig_ex*v;

        double nu_ion =
            nn*sig_ion*v;

        double nuReal =
            nu_el +
            nu_ex +
            nu_ion;

        //----------------------------------------------
        // null collision
        //----------------------------------------------

        if(rnd() > nuReal/nuMax)
            continue;

        //----------------------------------------------
        // choose collision type
        //----------------------------------------------

        double r =
            rnd()*nuReal;

        //--------------------------------------------------
        // elastic
        //--------------------------------------------------

        if(r < nu_el)
        {
            ep.vel =
                randomScatter2(
                    ep.vel,
                    300.0,
                    neutrals.mass
                );
        }

        //--------------------------------------------------
        // excitation
        //--------------------------------------------------

        else if(r < nu_el+nu_ex)
        {
            const double Eex = 8.0;

            if(Ee > Eex)
            {
                double Enew =
                    Ee-Eex;

                double scale =
                    sqrt(Enew/Ee);

                ep.vel *= scale;
            }
        }

        //--------------------------------------------------
        // ionization
        //--------------------------------------------------

        else
        {
            const double Eion = 15.8;

            if(Ee <= Eion)
                continue;

            //------------------------------------------
            // primary electron
            //------------------------------------------

            double Eremain =
                Ee-Eion;

            double Eprimary =
                0.5*Eremain;

            double scale =
                sqrt(Eprimary/Ee);

            ep.vel *= scale;

            //------------------------------------------
            // secondary electron
            //------------------------------------------

            Particle sec;

            sec.pos =
                ep.pos;

            sec.vel =
                sampleIonizationElectron();

            sec.STATUS =
                PARTICLE_ACTIVE;

            newElectrons.push_back(sec);

            //------------------------------------------
            // ion source
            //------------------------------------------

            int i =
                floor(lc[0]);

            int k =
                floor(lc[2]);

            int kg =
                k+1;

            double r1 =
                std::max(
                    0.0,
                    (i-0.5)*domain.dh[0]
                );

            double r2 =
                (i+0.5)
                *
                domain.dh[0];

            double cellVol =
                Const::PI
                *
                (r2*r2-r1*r1)
                *
                domain.dh[2];

            ions.source(i,kg)
                +=
                electrons.spwt
                /
                (
                    cellVol
                );
        }
    }

    //--------------------------------------------------
    // append new electrons
    //--------------------------------------------------

    for(auto& p : newElectrons)
    {
        electrons.particles.push_back(p);
    }
}


void MCC_CEX ::	apply(double dt)
{
    int Nr = domain.Nr;

    int Nz = domain.Nz_phys;

    double sigma =
        1e-18;

    for(int k=0;k<Nz;k++)
    {
        int kg = k+1;

        for(int i=0;i<Nr;i++)
        {
            double nn =
                neutrals.den(i,kg);

            double ur =
                ions.ur(i,kg);

            double uz =
                ions.uz(i,kg);

            double ui =
                sqrt(
                    ur*ur +
                    uz*uz
                );

            double nu =
                nn*sigma*ui;

            double factor =
                exp(-nu*dt);

            ions.ur(i,kg) *= factor;

            ions.uz(i,kg) *= factor;
        }
    }
}



double DSMC_MEX :: evalSigma(double g_rel)
{
    return
        Const::PI* c[0]*c[0] *
        pow(
            c[2]/(g_rel*g_rel),
            c[1]-0.5
        	)
        /
        c[3];
}

/* collides two particles*/
void DSMC_MEX :: collide(
					vector3d& vel1,	vector3d& vel2,
					double mass1, 	double mass2
    	)
{
    //----------------------------------
    // center-of-mass velocity
    //----------------------------------

    vector3d cm =
        (vel1*mass1 + vel2*mass2)
        *
        (1.0/(mass1+mass2));

    //----------------------------------
    // relative velocity
    //----------------------------------

    vector3d cr =
        vel1 - vel2;

    double cr_mag =
        cr.norm();

    //----------------------------------
    // isotropic scattering
    //----------------------------------

    double cos_chi =
        2.0*rnd() - 1.0;

    double sin_chi =
        sqrt(
            1.0 -
            cos_chi*cos_chi
        );

    double phi =
        2.0*
        Const::PI*
        rnd();

    //----------------------------------
    // new relative velocity
    //----------------------------------

    vector3d cr_new;

    cr_new.x =
        cr_mag*cos_chi;

    cr_new.y =
        cr_mag*sin_chi*cos(phi);

    cr_new.z =
        cr_mag*sin_chi*sin(phi);

    //----------------------------------
    // post-collision velocities
    //----------------------------------

    vel1 =
        cm
        +
        cr_new*
        (mass2/(mass1+mass2));

    vel2 =
        cm
        -
        cr_new*
        (mass1/(mass1+mass2));
}


/*DSMC*/
void DSMC_MEX::apply(double dt)
{
    double lc[3];

    const int Nr =
        domain.Nr;

    const int Nz =
        domain.Nz_phys;

    const int n_cells =
        Nr * Nz;

    std::vector<Particle*>* parts_in_cell =
        new std::vector<Particle*>[n_cells];

    //--------------------------------------------------
    // sort particles
    //--------------------------------------------------

    for (Particle& part : species.particles)
    {
        if (
            part.STATUS
            != PARTICLE_ACTIVE
        )
            continue;

        domain.Posv3d2lc(
            lc,
            part.pos
        );

        int i =
            (int)floor(lc[0]);

        int k =
            (int)floor(lc[2]);

        if (
            i < 0 ||	i >= Nr ||
            k < 0 || k >= Nz
        	)
            continue;

        int c =
            i + k*Nr;

        parts_in_cell[c]
            .push_back(&part);
    }

    //--------------------------------------------------
    // DSMC collisions
    //--------------------------------------------------

    double sigma_cr_max_temp =
        sigma_cr_max;

    int num_cols = 0;

    double Fn =
        species.spwt;

    for (int k=0;k<Nz;k++)
    {
        for (int i=0;i<Nr;i++)
        {
            int c =
                i + k*Nr;

            auto& parts =
                parts_in_cell[c];

            int np =
                (int)parts.size();

            if (np < 2)
                continue;

            //------------------------------------------
            // cylindrical cell volume
            //------------------------------------------

            double r =
                std::max(
                    (i + 0.5)
                    *
                    domain.dh[0],
                    0.5*
                    domain.dh[0]
                );

            double dV =
                2.0*
                Const::PI*
                r*
                domain.dh[0]*
                domain.dh[2];

            //------------------------------------------
            // Bird NTC
            //------------------------------------------

            double ng_f =
                0.5 *
                np *
                np *
                Fn *
                sigma_cr_max *
                dt /
                dV;

            int ng =
                (int)(ng_f + 0.5);

            for (int g=0;
                 g<ng;
                 g++)
            {
                int p1 =
                    (int)(
                        rnd()*np
                    );

                int p2;

                do
                {
                    p2 =
                        (int)(
                            rnd()*np
                        );
                }
                while (p2 == p1);

                vector3d cr_vec =
                    parts[p1]->vel
                    -
                    parts[p2]->vel;

                double cr =
                    cr_vec.norm();

                double sigma =
                    evalSigma(cr);

                double sigma_cr =
                    sigma * cr;

                sigma_cr_max_temp =
                    std::max(
                        sigma_cr_max_temp,
                        sigma_cr
                    );

                double P =
                    sigma_cr
                    /
                    std::max(
                        sigma_cr_max,
                        1e-30
                    );

                if (P > rnd())
                {
                    collide(
                        parts[p1]->vel,
                        parts[p2]->vel,
                        species.mass,
                        species.mass
                    );

                    num_cols++;
                }
            }
        }
    }

    delete[] parts_in_cell;

    if (num_cols > 0)
    {
        sigma_cr_max =
            sigma_cr_max_temp;
    }
}



vector3d randomScatter(
    const vector3d& vin
)
{
    double vmag = vin.norm();

    //------------------------------------
    // isotropic direction
    //------------------------------------

    double mu =
        2.0*rnd() - 1.0;

    double phi =
        2.0*Const::PI*rnd();

    double sint =
        sqrt(1.0 - mu*mu);

    vector3d vout;

    vout.x =
        vmag * sint * cos(phi);

    vout.y =
        vmag * sint * sin(phi);

    vout.z =
        vmag * mu;

    return vout;
}

vector3d randomScatter2(
    const vector3d& ve,	double Tn,	double mn
	)
{
    //------------------------------------
    // sample neutral velocity
    //------------------------------------
    vector3d vn =
        sampleMaxwellianVelocity(
            Tn,	mn
        );
    //------------------------------------
    // relative velocity
    //------------------------------------

    vector3d g =
        ve - vn;

    double gmag =
        g.norm();

    //------------------------------------
    // isotropic scattering
    //------------------------------------

    double mu =
        2.0*rnd() - 1.0;

    double phi =
        2.0*Const::PI*rnd();

    double sint =
        sqrt(1.0 - mu*mu);

    vector3d gnew;

    gnew.x =	gmag*sint*cos(phi);

    gnew.y =	gmag*sint*sin(phi);

    gnew.z =	gmag*mu;
    //------------------------------------
    // back to lab frame
    //------------------------------------
    return vn + gnew;
}