#ifndef CURRENTLOOP_H
#define CURRENTLOOP_H

#include <cmath>
#include "Utils.h"

class CurrentLoop
{
public:

    CurrentLoop();
    ~CurrentLoop() = default;

    //----------------------------------
    // Electromagnetic fields
    //----------------------------------

    vector3d solveB(const vector3d& position) const;

    vector3d solveA(const vector3d& position) const;

    //----------------------------------
    // Coil parameters
    //----------------------------------

    vector3d center;      // [m]

    vector3d normal;      // unit vector

    double radius;        // [m]

    double current;       // [A]

    //----------------------------------
    // Numerical integration
    //----------------------------------

    int integrationSteps;
	
private:
	
    void buildLocalBasis(
        vector3d& e1,
        vector3d& e2
    ) const;
};

#endif