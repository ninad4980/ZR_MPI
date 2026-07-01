#ifndef RFFIELDSOLVER_H
#define RFFIELDSOLVER_H

#include "Domain.h"
#include "RFCoilAssembly.h"

class RFFieldSolver
{
public:

    RFFieldSolver(
        Domain& domain,
        RFCoilAssembly& antenna
    );

    void computeA();

    void computeB();

    void computeERF();

//	void update(double time);
	
private:

    Domain& domain;

    RFCoilAssembly& antenna;
};

#endif