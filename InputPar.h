#ifndef INPUTPAR_H
#define INPUTPAR_H

#include <string>
#include <map>

class InputPar
{
public:

    //--------------------------------------------------
    // geometry
    //--------------------------------------------------

	 struct Domain_t
    {
        int Nr;
        int Nzee;
        int Ntheta;

        double Radius;
        double Length;
        
        double posFC;
    } 
    	Domain;
    
    struct Chamber_t
    {
    double WallThick;

    double InletAperture;
    double InletRadius;
    double InletStart;
    double InletLen;

    double Radius;
    double Length;
    double ChamberInletRadius;
    double ChamberExtractRadius;

	double phiChamber;
    double phiExtract;

    double ExtractGap;
    double ExtractThick;
    double ExtractrAperture;
    double ExtractRout;

    double InsulatorThick;
    double InsulatorRout;
    double InsulatorAperture;
     //--------------------------------------------------
    // derived geometry
    //--------------------------------------------------

    double ChamberStart;
    double ChamberEnd;

    double ExtractPos;

    double InsulatorPos;

    double inletRadius;
    double zInject;

    } 
    	Chamber;
    	


    

    //--------------------------------------------------
    // simulation
    //--------------------------------------------------
	 struct Plasma_t
    {
    double initial_E_Density;
    double initialIonDensity;

    double kTe;
    double kTi_O2;
	}
	Plasma;
	
	struct Gas_t
	{
    double inletPressure_N2;
    double inletPressure_O2;

    double inletTemperature;

    double bufferDensity;
	}
	Gas;
		

    //--------------------------------------------------
    // methods
    //--------------------------------------------------
	struct RF_t
    {
        bool Enabled;

        double Power;
        double Frequency;
        double Current;
		double Phase;
		
        int Turns;
        int Layers;
		double LayerSpacing;
		
        double Radius;
        double Pitch;
        double StartZ;
    } RF;

    struct StaticField_t
    {
        bool Enabled;

        double Current;

        int Turns;

        double Radius;

        double Spacing;

        double TargetB;

        int Configuration;
    } StaticField;

    struct Numerics_t
    {
        double dtElectron;

        double dtIon;

        double dtNeutral;
    } Numerics;

    InputPar();

    bool readInput(
        const std::string &filename
    );

private:

    std::map<std::string,double> values;

    double get(
        const std::string &key
    ) const;
};

extern InputPar inputPar;

#endif