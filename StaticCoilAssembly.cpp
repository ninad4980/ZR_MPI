#include "StaticCoilAssembly.h"

using namespace Const;

///////////////////////////////////////////////////////////
//
// Constructor
//
///////////////////////////////////////////////////////////

StaticCoilAssembly::StaticCoilAssembly(const InputPar& input)
{
    center = {0.0,0.0,0.0};

    normal = {0.0,0.0,1.0};

	current = input.StaticField.Current;

    turns = input.StaticField.Turns;

    coilRadius =  input.StaticField.Radius;

    coilSpacing = input.StaticField.Spacing;

    targetB = input.StaticField.TargetB;
    
    chamberLength =	input.Chamber.Length;
    
    chamberRadius = input.Chamber.Radius;
    
    createHelmholtz();
    
    switch(input.StaticField.Configuration)
		{
		case 0:

	    createSingleCoil();
    	break;

		case 1:

    	createHelmholtz();
	    break;

		case 2:

    	createMirror();
    	break;

		default:

    	createHelmholtz();
		}
    
    
}

///////////////////////////////////////////////////////////
//
// Single coil
//
///////////////////////////////////////////////////////////

void StaticCoilAssembly::createSingleCoil()
{
    coils.clear();

    CurrentLoop loop;

    loop.center = center;

    loop.normal = normal;

    loop.radius = coilRadius;

    loop.current = current*turns;

    coils.push_back(loop);
}

///////////////////////////////////////////////////////////
//
// Helmholtz coils
//
///////////////////////////////////////////////////////////

void StaticCoilAssembly::createHelmholtz()
{
    coils.clear();

    CurrentLoop c1;
    CurrentLoop c2;

    c1.radius = coilRadius;
    c2.radius = coilRadius;

    c1.current = current*turns;
    c2.current = current*turns;

    c1.normal = normal;
    c2.normal = normal;

    c1.center = center;
    c2.center = center;

    c1.center.z -= 0.5*coilSpacing;

    c2.center.z += 0.5*coilSpacing;

    coils.push_back(c1);
    coils.push_back(c2);
}

///////////////////////////////////////////////////////////
//
// Magnetic mirror
//
///////////////////////////////////////////////////////////

void StaticCoilAssembly::createMirror()
{
    coils.clear();

    CurrentLoop c1;
    CurrentLoop c2;

    c1.radius = coilRadius;
    c2.radius = coilRadius;

    c1.current = current*turns;
    c2.current = current*turns;

    c1.normal = normal;
    c2.normal = normal;

    c1.center = center;
    c2.center = center;

    c1.center.z = center.z - 0.45*chamberLength;

    c2.center.z = center.z + 0.45*chamberLength;

    coils.push_back(c1);
    coils.push_back(c2);
}

///////////////////////////////////////////////////////////
//
// Total magnetic field
//
///////////////////////////////////////////////////////////

vector3d StaticCoilAssembly::solveB(
    const vector3d& position
) const
{
    vector3d B;

    B = {0.0,0.0,0.0};

    for(const auto& coil : coils)
    {
        B += coil.solveB(position);
    }

    return B;
}