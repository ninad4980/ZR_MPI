#ifndef RFCOILASSEMBLY_H
#define RFCOILASSEMBLY_H

#include <vector>

#include "Utils.h"
#include "InputPar.h"
#include "CurrentLoop.h"

class RFCoilAssembly
{
public:

    RFCoilAssembly(const InputPar& input);
	
    ~RFCoilAssembly() = default;

    //------------------------------------
    // Geometry
    //------------------------------------

    void build();

    //------------------------------------
    // Electromagnetic fields
    //------------------------------------

    vector3d solveB(
        const vector3d& pos
    ) const;

    vector3d solveA(
        const vector3d& pos
    ) const;

    //------------------------------------
    // Coil parameters
    //------------------------------------

    std::vector<CurrentLoop> loops;

    vector3d center;

    vector3d normal;

    //------------------------------------
    // User inputs
    //------------------------------------

    int turns;

	int layers;
	
    double radius;

    double pitch;

    double current;

    double power;

    double frequency;

    double phase;

	
};

#endif