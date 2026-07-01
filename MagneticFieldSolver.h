#ifndef MAGNETICFIELDSOLVER_H
#define MAGNETICFIELDSOLVER_H

#include "Domain.h"
#include "StaticCoilAssembly.h"

class MagneticFieldSolver
{
public:

    MagneticFieldSolver(
        Domain& domain,
        StaticCoilAssembly& coils
    );

    //----------------------------------
    // compute static magnetic field
    //----------------------------------

    void compute();

private:

    Domain& domain;

    StaticCoilAssembly& coils;
};

#endif