#include "InputPar.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <algorithm>
InputPar inputPar;

InputPar::InputPar()
{
    
}
double InputPar::get(
    const std::string &key
) const
{
    auto it = values.find(key);

    if(it == values.end())
    {
        std::cerr
            << "Missing parameter: "
            << key
            << std::endl;

        std::exit(EXIT_FAILURE);
    }

    return it->second;
}

bool InputPar::readInput(
    const std::string &filename
)
{
    std::ifstream in(filename);

    if(!in)
    {
        std::cerr
            << "Cannot open "
            << filename
            << std::endl;

        return false;
    }

    std::string line;

    while(std::getline(in,line))
    {
        if(line.empty())
            continue;

        if(line[0]=='#')
            continue;
		
		std::replace(
        line.begin(),
        line.end(),
        '=',
        ' '
    	);
    	
        std::stringstream ss(line);

        std::string key;
        double value;

        ss >> key >> value;

        if(ss)
        {
            values[key] = value;
        }
    }
    //--------------------------------------------------
    // Domain
    //--------------------------------------------------
	Domain.Nr          = get("Nr");
	Domain.Ntheta          = get("Ntheta");
	Domain.Nzee          = get("Nzee");
	Domain.Radius          = get("DomainCylR");
	Domain.Length         = get("DomainCylLen");
	Domain.posFC              = get("posFC");
    //--------------------------------------------------
    // geometry
    //--------------------------------------------------

    Chamber.WallThick          = get("WallThick");

    Chamber.InletAperture      = get("InletAperture");
    Chamber.InletRadius        = get("InletRadius");
    Chamber.InletStart         = get("ChamberInletStart");
    Chamber.InletLen           = get("InletLen");

    Chamber.Radius      = get("ChamberRadius");
    Chamber.Length         = get("ChamberLen");

    Chamber.ChamberExtractRadius =
        get("ChamberExtractRadius");

    Chamber.phiChamber         = get("phiChamber");
    Chamber.phiExtract         = get("phiExtract");

    Chamber.ExtractGap         = get("ExtractGap");

    Chamber.ExtractThick       = get("ExtractThick");

    Chamber.ExtractrAperture   =
        get("ExtractrAperture");

    

    //--------------------------------------------------
    // simulation
    //--------------------------------------------------

    Plasma.initial_E_Density =
        get("initial_E_Density");

    Plasma.initialIonDensity =
        get("initialIonDensity");

    Plasma.kTe =
        get("kTe");

    Plasma.kTi_O2 =
        get("kTi_O2");

    Gas.inletPressure_O2 =
        get("inletPressure_O2");

    Gas.inletPressure_N2 =
        get("inletPressure_N2");

    Gas.inletTemperature =
        get("inletTemperature");

    Gas.bufferDensity =
        get("bufferDensity");

    //--------------------------------------------------
    // derived values
    //--------------------------------------------------

    Chamber.ChamberInletRadius =
			Chamber.InletRadius;

    Chamber.ChamberStart =
   		 Chamber.InletStart
    	+ Chamber.WallThick
    	+ Chamber.InletLen;

    Chamber.ChamberEnd =
    	Chamber.ChamberStart
    	+ Chamber.WallThick
    	+ Chamber.Length
    	+ Chamber.WallThick;

    Chamber.ExtractPos =
    	Chamber.ChamberEnd
    	+ Chamber.ExtractGap;

   Chamber.ExtractRout =
    	Chamber.Radius
    	+ Chamber.WallThick;

    Chamber.InsulatorPos =
   		 Chamber.ChamberEnd;

    Chamber.InsulatorThick =
   		 Chamber.ExtractGap;

   Chamber.InsulatorRout =
    	Chamber.Radius
    	+ Chamber.WallThick;

    Chamber.InsulatorAperture =
   		 Chamber.InsulatorRout
    		- 4e-3;

    Chamber.inletRadius =
   		 Chamber.InletRadius
    		-1e-3;

    Chamber.zInject =
  	  Chamber.InletStart
   		 +4e-3;
	//--------------------------------------------------
	// RF
	//--------------------------------------------------

	RF.Enabled    = static_cast<bool>(get("RFEnabled"));

	RF.Power      = get("rfPower");

	RF.Frequency  = get("rfFrequency");

	RF.Current    = get("rfCurrent");

	RF.Turns      = static_cast<int>(get("rfTurns"));

	RF.Layers     = static_cast<int>(get("rfLayers"));

	RF.Radius     = get("rfRadius");

	RF.Pitch      = get("rfPitch");

	RF.StartZ     = get("rfStartZ");
	
	RF.LayerSpacing = get("rfLayerSpacing");
	//--------------------------------------------------
	// Static field
	//--------------------------------------------------

	StaticField.Enabled =
 		   static_cast<bool>(
  	      get("StaticFieldEnabled")
  	 	 );

	StaticField.Current =
 		   get("StaticCurrent");

	StaticField.Turns =
 		   static_cast<int>(
 	       get("StaticTurns")
    	);

	StaticField.Radius =
   		 get("StaticRadius");

	StaticField.Spacing =
    	get("StaticSpacing");

	StaticField.TargetB =
  		  get("TargetB");

	StaticField.Configuration =
  		  static_cast<int>(
     	   get("StaticConfiguration")
   		 );
   		 
   		 
   	Numerics.dtElectron =
   		 get("dtElectron");

	Numerics.dtIon =
  		  get("dtIon");

	Numerics.dtNeutral =
    		get("dtNeutral");	 
    		
    		
    return true;
}