#ifndef STATICCOILASSEMBLY_H
#define STATICCOILASSEMBLY_H

#include <vector>

#include "CurrentLoop.h"
#include "Utils.h"
#include "InputPar.h"
class StaticCoilAssembly
{
public:

    StaticCoilAssembly(const InputPar& input);
    ~StaticCoilAssembly() = default;

    //---------------------------------------
    // Build common coil configurations
    //---------------------------------------

    void createSingleCoil();

    void createHelmholtz();

    void createMirror();

    //---------------------------------------
    // Field evaluation
    //---------------------------------------

    vector3d solveB(const vector3d& position) const;

    //---------------------------------------
    // Coil list
    //---------------------------------------

    std::vector<CurrentLoop> coils;

    //---------------------------------------
    // User parameters
    //---------------------------------------

    vector3d center;

    vector3d normal;

    double chamberRadius;

    double chamberLength;

    double coilRadius;

    double current;

    int turns;

    double coilSpacing;
    
    double targetB;

};

#endif