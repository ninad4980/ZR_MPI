#ifndef _DOMAINBUILDER_H
#define _DOMAINBUILDER_H

#include "Domain.h"
#include "InputPar.h"

class DomainBuilder
{
public:

    static Domain create(
        const InputPar& input,
        int mpi_rank,
        int mpi_size
    );
};


#endif