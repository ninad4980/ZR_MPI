#include "GeometryBuilder.h"

#include "Geometry.h"

void GeometryBuilder::build(
    Domain& domain,
    const InputPar& input
)
{
    GeometryManual geo(domain);

    //--------------------------------------------------
    // Plasma chamber
    //--------------------------------------------------

    geo.plasmaChamber();

    //--------------------------------------------------
    // Gas inlet
    //--------------------------------------------------

    geo.AddObject_Inlet();

    //--------------------------------------------------
    // Extraction electrode
    //--------------------------------------------------

    geo.AddObject_electrode(
        inputPar.Chamber.ExtractPos,
        inputPar.Chamber.ExtractThick,
        inputPar.Chamber.ExtractrAperture,
        inputPar.Chamber.ExtractRout,
        inputPar.Chamber.phiExtract,
        OBJ_ELECTRODE
    );

    //--------------------------------------------------
    // Insulator
    //--------------------------------------------------

    geo.AddObject_electrode(
        inputPar.Chamber.InsulatorPos,
        inputPar.Chamber.InsulatorThick,
        inputPar.Chamber.InsulatorAperture,
        inputPar.Chamber.InsulatorRout,
        0.0,
        OBJ_DIELECTRIC
    );
}