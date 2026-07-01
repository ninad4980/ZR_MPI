#include "DomainBuilder.h"

#include "Domain.h"
#include "InputPar.h"

Domain DomainBuilder::create(
    const InputPar& input,
    int mpi_rank,
    int mpi_size
)
{
    Domain domain(
        input.Domain.Nr,
        input.Domain.Ntheta,
        input.Domain.Nzee,
        (input.Domain.Nzee - 1)/mpi_size + 3
    );

    const double dr =
        input.Domain.Radius /
        (input.Domain.Nr - 1);

    const double dz =
        input.Domain.Length /
        (input.Domain.Nzee - 1);

    domain.setOrigin(
        0.0,
        0.0,
        0.0
    );

    domain.setSpacing(
        dr,
        1.0,
        dz
    );

    domain.initMPIDomain(
        mpi_rank,
        mpi_size
    );

    return domain;
}