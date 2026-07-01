#include "Geometry.h"
#include "InputPar.h"
/*****************************************************************/
/**							    	**/
	/** Following code is for manually setting geometry **/
/**								**/
/*****************************************************************/

/** Only defines region **/
void GeometryManual::plasmaChamber()
{
    Field &phi = domain.phi;

    for (int k=0; k<Nz_phys; k++)
    {
        for (int i=0; i<Nr; i++)
        {
            double r[3];

            domain.pos(r,i,0,k);

            r[1] = 0.0;

            int id = i + k*Nr;

            //--------------------------------------------------
            // default
            //--------------------------------------------------

            domain.object[id] = OBJ_VOID;

            //--------------------------------------------------
            // BACK WALL
            //--------------------------------------------------

            if (
                r[2] >= inputPar.Chamber.ChamberStart &&
                r[2] <= inputPar.Chamber.ChamberStart + 
                		inputPar.Chamber.WallThick
            	)
            {
                if (r[0] < inputPar.Chamber.ChamberInletRadius)
                {
                    domain.object[id] =
                        OBJ_PLASMA;
                }
                else if (
                    r[0] >= inputPar.Chamber.ChamberInletRadius &&
                    r[0] <= inputPar.Chamber.Radius + 
                    		inputPar.Chamber.WallThick
               		)
                {
                    domain.object[id] =
                        OBJ_WALL;

                    phi(i,k+1) = inputPar.Chamber.phiChamber;
                }
            }

            //--------------------------------------------------
            // CYLINDER
            //--------------------------------------------------

            else if (
                r[2] > inputPar.Chamber.ChamberStart + inputPar.Chamber.WallThick &&
                r[2] < inputPar.Chamber.ChamberStart + inputPar.Chamber.WallThick
                		 + inputPar.Chamber.Length
           		)
            {
                if (r[0] < inputPar.Chamber.Radius)
                {
                    domain.object[id] =
                        OBJ_PLASMA;
                }
                else if (
                    r[0] >= inputPar.Chamber.Radius &&
                    r[0] <= inputPar.Chamber.Radius 
                    		+ inputPar.Chamber.WallThick
                )
                {
                    domain.object[id] =
                        OBJ_WALL;

                    phi(i,k+1) = inputPar.Chamber.phiChamber;
                }
            }

            //--------------------------------------------------
            // FRONT WALL
            //--------------------------------------------------

            else if (
                r[2] >= inputPar.Chamber.ChamberStart + inputPar.Chamber.WallThick 
                		+ inputPar.Chamber.Length &&
                r[2] <  inputPar.Chamber.ChamberStart + inputPar.Chamber.WallThick +
                 	inputPar.Chamber.Length + inputPar.Chamber.WallThick
            	)
            {
                if (r[0] < inputPar.Chamber.ChamberExtractRadius )
                {
                    domain.object[id] =
                        OBJ_PLASMA;
                }
                else if(r[0]>= inputPar.Chamber.ChamberExtractRadius && 
                			r[0]<=inputPar.Chamber.ExtractRout )
                {
                    domain.object[id] =
                        OBJ_WALL;

                    phi(i,k+1) = inputPar.Chamber.phiChamber;
                }
            }

            //--------------------------------------------------
            // EXTRACTION REGION
            //--------------------------------------------------

            else if (
                r[2] >=
                inputPar.Chamber.ChamberEnd
            	)
            {
                domain.object[id] =
            	        OBJ_EXTRACTION;
            }
        // end of i,k loop  
        }
    }

// end of function 
}
/** Cylindrical object.data : starting position, length, inner radius, outer radius
potential on cylinder, object.data index **/

void GeometryManual::AddObject_Inlet()
{
    Field &phi = domain.phi;

    for (int k=0; k<Nz_phys; k++)
    {
        for (int i=0; i<Nr; i++)
        {
            //--------------------------------------------------
            // node position
            //--------------------------------------------------

            double r[3];

            domain.pos(r,i,0,k);

            r[1] = 0.0;

            int id = i + k*Nr;

			//--------------------------------------------------
			// inlet back plate with aperture
			//--------------------------------------------------

			if (
			    r[2] <= inputPar.Chamber.WallThick
				)
				{
		    	if (r[0] < inputPar.Chamber.InletAperture)
    				{
			        domain.object[id] =
            				OBJ_PLASMA;
    				}
		     	if (
        			r[0] >= inputPar.Chamber.InletAperture && 
        					r[0]<= (inputPar.Chamber.InletRadius+inputPar.Chamber.WallThick)
    				)
    				{
			        domain.object[id] =
        		    OBJ_WALL;

			        phi(i,k+1) =
        			    inputPar.Chamber.phiChamber;
    				}
				}
            //--------------------------------------------------
            // inlet axial region
            //--------------------------------------------------

            if (
				 r[2] >  (inputPar.Chamber.WallThick)	
                &&
                r[2] <= (inputPar.Chamber.ChamberStart)
           			)
            	{
                //--------------------------------------------------
                // inlet opening
                //--------------------------------------------------

                if (r[0] < inputPar.Chamber.InletRadius)
                {
                    domain.object[id] =
                        OBJ_PLASMA;
                }

                //--------------------------------------------------
                // inlet wall
                //--------------------------------------------------

                else if (
                    r[0] >= inputPar.Chamber.InletRadius
                    &&
                    r[0] <= inputPar.Chamber.InletRadius + inputPar.Chamber.WallThick
                )
                {
                    domain.object[id] =
                        OBJ_WALL;

                    phi(i,k+1) =
                        inputPar.Chamber.phiChamber;
                }
                
                
            }
        }
    }
}

/** Disc like object.data with a hole : starting Position, thickness, radius of hole
outer radius of disc , Potential on electrode, object.data index  **/
void GeometryManual::AddObject_electrode(
    double posE,
    double dE,
    double rAperture,
    double rOut,
    double phi_electrode,
    ObjectType type
)
{
    Field &phi = domain.phi;

    for (int k=0; k<Nz_phys; k++)
    {
        for (int i=0; i<Nr; i++)
        {
            //--------------------------------------------------
            // node position
            //--------------------------------------------------

            double r[3];

            domain.pos(r,i,0,k);

            r[1] = 0.0;

            int id = i + k*Nr;

            //--------------------------------------------------
            // electrode axial extent
            //--------------------------------------------------

            if (
                r[2] >= posE
                &&
                r[2] <  posE + dE
            )
            {
                //--------------------------------------------------
                // annular conducting region
                //--------------------------------------------------

                if (
                    r[0] >= rAperture
                    &&
                    r[0] <  rOut
                )
                {
                    domain.object[id] =
                        type;

                    phi(i,k+1) =
                        phi_electrode;
                }
            }
        }
    }
}
