#ifndef _GEOMETRYBUILDER_H
#define _GEOMETRYBUILDER_H

#include "Domain.h"
#include "InputPar.h"

class GeometryBuilder
{
public:

    static void build(
    Domain& domain,
    const InputPar& input
	);
    
};

#endif