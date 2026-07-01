/***************************************************************************
 *   Copyright (C) 2009 by Ninad Joshi   *
 ***************************************************************************/
/*Field is a container for mesh node data
It contains functions for scatter/gather, and division by volume
definition of vector3d function
*/


#include <math.h>
#include <iostream>
#include <fstream>
#include <mpi.h>
#include <sstream>
#include <iomanip>
#include <random>
#include <chrono>
#include <thread>
#include <string>
#include <cstring>
#include <omp.h>
#include<sys/stat.h>
#include <memory>
#include "Core.h"
#include "GeometryModule.h"
#include "PhysicsModule.h"
#include "EMModule.h"
#include "IOModule.h"

/*function prototypes*/

using namespace Const;		//to avoid having to write Const::ME

/*simulation parameters*/


/*globals*/
int mpi_rank=0;
int nthreads =  omp_get_num_threads(), thread_num = omp_get_thread_num();

/*helper function to print message only on root and abort*/
void error(std ::string msg)
{
	if (mpi_rank==0) std :: cerr<<msg<<std ::endl;
	exit(-1);
}





/*main function*/
int main(int n_args, char *args[])
{

	/* Initialize MPI */
	MPI_Init(NULL, NULL);

	/*get number of processes from MPI*/
	int mpi_size;
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

	/*get our rank*/
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

 	/*initialize domain with global parameters*/

	std :: cout<<"rank = "<< mpi_rank<<"\n";
	if (mpi_rank==0)
		{	if(mkdir("results",0777)==-1)
				{std::cerr<<"Msg: "<<strerror(errno)<<std::endl;}
			else 
				{std::cout<<"Directory created \n";}	
		}
	MPI_Barrier( MPI_COMM_WORLD);			
	omp_set_dynamic(0);
	omp_set_num_threads(2);
	
	  if(!inputPar.readInput("input.dat"))
        return -1;

    std::cout
        << "Input file loaded."
        << std::endl;

	/*capture starting time*/
	//auto clock_start = chrono::high_resolution_clock::now();
	/** The main grid object where fields are calculated (r,theta,zee) **/
	/* Domain .h .cpp				              */
	/* Initilize with ni,nj,nk, nklocal , in meters		      */
	/* Only decomposed in nk/z direction 		              */
	/*****************************************************************/

Domain domain =
    DomainBuilder::create(
        inputPar,
        mpi_rank,
        mpi_size
    );

	GeometryBuilder::build(
    domain,
    inputPar
	);
					
	MPI_Barrier( MPI_COMM_WORLD);
	
	StaticCoilAssembly staticCoils(inputPar);

	MagneticFieldSolver magneticSolver(
		    domain,
    		staticCoils
			);

	magneticSolver.compute();

	MPI_Barrier( MPI_COMM_WORLD);

	RFCoilAssembly rfCoils(inputPar);

	RFFieldSolver rfSolver(
    		domain,
		    rfCoils
			);

	rfSolver.computeA();
	rfSolver.computeB();
	rfSolver.computeERF();
	
	MPI_Barrier( MPI_COMM_WORLD);
	std::string str;
	char *cstr = NULL;

	Output output(domain); 

	output.write2file_Object();

	MPI_Barrier( MPI_COMM_WORLD);
		
	//--------------------------------------------------
	// set up particle species
	//--------------------------------------------------

	std::vector<std::unique_ptr<Species>> speciesList;

	speciesList.reserve(10);
	//--------------------------------------------------
	// Neutral DSMC species
	//--------------------------------------------------

	auto O2 =
    	std::make_unique<KineticSpecies>(
    	    domain,
    	    NEUTRAL,
    	    "O2",
    	    32.0 * AMU,
       		0.0,
       		1e9
    	);

	KineticSpecies* O2_ptr = O2.get();

	speciesList.push_back(std::move(O2));

	//--------------------------------------------------
	// Electron PIC species
	//--------------------------------------------------

	auto electrons =
    	std::make_unique<KineticSpecies>(
    	    domain,
    	    ELECTRON,
    	    "e-",
    	    ME,
    	    -QE,
    	    1e6
    	);

	KineticSpecies* e_ptr = electrons.get();
	
	speciesList.push_back(std::move(electrons));
	
	//--------------------------------------------------
	// Ions fluid species
	//--------------------------------------------------

	auto ions =
    std::make_unique<FluidSpecies>(
        domain,
        ION,
        "O2+",
        32.0 * AMU,
        QE
    	);
    
	FluidSpecies* O2plus_ptr = ions.get();
	
	speciesList.push_back(std::move(ions));
	//--------------------------------------------------
	// Set collisions
	//--------------------------------------------------
	CollisionManager collisions;


	collisions.interactions.push_back(
    			std::make_unique<DSMC_MEX>(
        			*O2_ptr,
        				domain
    			)
			);
	collisions.interactions.push_back(
    			std::make_unique<ElectronNeutralMCC>(
        		*e_ptr,
			        *O2_ptr,
				        *O2plus_ptr,
        					domain
    			)
			);
	collisions.interactions.push_back(
    			std::make_unique<MCC_CEX>(
        		*O2plus_ptr,
        			*O2_ptr,
        					domain
    			)
			);
	//--------------------------------------------------
	// Set small time step for each species
	//--------------------------------------------------
	e_ptr->dt      = inputPar.Numerics.dtElectron;

	O2plus_ptr->dt = inputPar.Numerics.dtIon;

	O2_ptr->dt     = inputPar.Numerics.dtNeutral;
	//--------------------------------------------------
	// Sources
	//--------------------------------------------------

	Source neutralSource(
	    domain,
	    *O2_ptr
		);

	Source eSource(
	    domain,
	    *e_ptr
		);
			
	//--------------------------------------------------
	// Initial particles
	//--------------------------------------------------

	neutralSource.injectNeutralInlet(
	    inputPar.Gas.inletPressure_O2,
	    inputPar.Gas.inletTemperature,
	    inputPar.Chamber.inletRadius-0.5e-03,
	    inputPar.Chamber.zInject,
	    O2_ptr->dt
		);

	eSource.samplePlasmaChamber_ConstDen(
	    e_ptr->spwt,
	    inputPar.Plasma.initial_E_Density,
	    inputPar.Plasma.kTe
		);

	initiateIonsPlasmaChamber_ConstDen(
  			  *O2plus_ptr,
    			domain,
			    inputPar.Plasma.initialIonDensity,
			    inputPar.Plasma.kTi_O2
			);
	//--------------------------------------------------
	// Movers
	//--------------------------------------------------

	 	DSMCMover neutralMover(
	 	   *O2_ptr,
	 	   domain
			);

		PICMover electronMover(
		    *e_ptr,
		    domain
			);
		
		FluidMover	O2plusMover(
			*O2plus_ptr,
			domain
			);
	//--------------------------------------------------
	// electrostatic solver
	//--------------------------------------------------
	PotentialSolver solver(
    	domain,
    	1000,
    	1e-4
		);
	
	//--------------------------------------------------
	// number densities
	//--------------------------------------------------

	for (auto &sp : speciesList)
	{
	    if (KineticSpecies *ks =
	        dynamic_cast<KineticSpecies*>(sp.get()))
	    {
	        ks->computeGasProperties();
	    }
	}

	//--------------------------------------------------
	// charge density
	//--------------------------------------------------

	domain.rho.clear();

	//--------------------------------------------------
	// timing
	//--------------------------------------------------

	double time = 0.0;

	double omega =
    	2.0*Const::PI*
		    inputPar.RF.Frequency;
    
	const int electronSubcycles =
    	static_cast<int>(
    	    O2_ptr->dt / e_ptr->dt
    	);
	const int ionSubcycles      = 
		static_cast<int>(
    	    O2_ptr->dt / O2plus_ptr->dt
    	);	
	//--------------------------------------------------
	// stages
	//--------------------------------------------------

		const int neutralFillSteps = 50001;
		const int plasmaSteps      = 2001;
	//--------------------------------------------------
	// STAGE 1
	// neutral filling only
	//--------------------------------------------------
	for(int neutralStep=0;
	    neutralStep<neutralFillSteps;
		    neutralStep++)
		{
    	time += O2_ptr->dt;
    	//------------------------------------------
    	// inject neutrals
    	//------------------------------------------
		neutralSource.injectNeutralInlet(
    		inputPar.Gas.inletPressure_O2,
        	inputPar.Gas.inletTemperature,
       		inputPar.Chamber.inletRadius-0.5e-03,
	        inputPar.Chamber.zInject,
	        O2_ptr->dt
			);
    	//-----------------------------------
    	// move neutrals
    	//-----------------------------------
    	neutralMover.move(0.0);
		//-----------------------------------
    		// DSMC collisions
    		//-----------------------------------
    	collisions.interactions[0]->apply(
      	  O2_ptr->dt
    		);   
		if (neutralStep % 1000 == 0)
    		{
    	//------------------------------------------
    	// diagnostics
    	//------------------------------------------
        	O2_ptr->computeGasProperties();

	        char fileName[128];

	        sprintf(
    	        fileName,
    	        "results/O2_fill_%05d.dat",
    	        neutralStep
    		    );
        	output.outputPart2(
        	    *O2_ptr,
        	    fileName
        		);
	        if (mpi_rank == 0)
		        {
    	        std::cout
        	        << "[Fill] step="
            	    << neutralStep
                	<< " neutrals="
                	<< O2_ptr->particles.size()
                	<< std::endl;
        		}
    		}
		}
	
	//--------------------------------------------------
	// STAGE 2
	// plasma evolution
	//--------------------------------------------------
	domain.rho.clear();

	solver.calcRhoC(speciesList);
	
	solver.solvePhi();

	solver.computeEF();
	
	for (int step=0;
		     step<plasmaSteps;
    		 step++)
		{
		time += O2_ptr->dt;
		
		double rfFactor =
		    std::cos(omega*time);
		
		O2plus_ptr->source.clear();

    	for(int esub=0; esub<electronSubcycles; esub++)
			{
		    double t =
    		    time
        		+ esub*e_ptr->dt;

    		double rfFactor =
        		std::cos(omega*t);

    		electronMover.move(rfFactor);

		    collisions.interactions[1]->apply(
        			e_ptr->dt);
			}	
    		
    	e_ptr->computeGasProperties();		
    	
    	//--------------------------------------
   		 // Ion fluid subcycling
    	 //--------------------------------------
    	 O2plus_ptr->source /=  O2plus_ptr->dt ;
	    for(int isub=0;
		        isub<ionSubcycles;
	   		    isub++)
	   		{
	   		double t =
			    time
			    + isub*O2plus_ptr->dt;

			double rfFactor =
			    std::cos(omega*t);

			O2plusMover.move(rfFactor);
	   		
	   		collisions.interactions[2]->apply(
		        O2plus_ptr->dt
   					 );
    		}
    		
    	O2_ptr->computeGasProperties();		
    	
    	//------------------------------------------
	    // rebuild rho
	    //------------------------------------------
    	solver.calcRhoC(speciesList);
    	//------------------------------------------
    	// solve electrostatic field
    	//------------------------------------------
    	solver.solvePhi();
    	
    	solver.computeEF();
	    //------------------------------------------
	    // diagnostics
	    //------------------------------------------
	    if (step % 20 == 0)
		    {
			char fileName[128];

        	//----------------------------------
        	// electrons
        	//----------------------------------

        	sprintf(
        	    fileName,
            	"results/e_%05d.dat",
            	step
        		);

	        output.outputPart2(
	            *e_ptr,
	            fileName
		        );
		
			
			sprintf(
    			fileName,
			    "results/ni_%05d.dat",
			    step
				);

			output.write2file_SField2d(
			    &O2plus_ptr->density,
			    fileName
				);
			sprintf(
  				fileName,
				"results/ne_%05d.dat",
			    step
				);

			output.write2file_SField2d(
			    &e_ptr->den,
			    fileName
				);
			sprintf(
			    fileName,
			    "results/nn_%05d.dat",
			    step
				);

			output.write2file_SField2d(
			    &O2_ptr->den,
			    fileName
				);
			
			sprintf(
			    fileName,
			    "results/rho_%05d.dat",
			    step
				);

			output.write2file_SField2d(
			    &domain.rho,
			    fileName
				);
			sprintf(
			    fileName,
			    "results/Phi_%05d.dat",
			    step
				);

			output.write2file_SField2d(
			    &domain.phi,
			    fileName
				);	
			if (mpi_rank == 0)
     		   {
            	std::cout
            	    << "[Plasma] step="
            	    << step
            	    << " neutrals="
            	    << O2_ptr->particles.size()
            	    << " electrons="
            	    << e_ptr->particles.size()
            	    << std::endl;
        		}		
			}
		}		
	if (mpi_rank==0)
    		{
		std :: cout<<"Running with "<<mpi_size<<" MPI domains"<<std ::endl;
    		}

	if (mpi_rank==0)
		{
		//cout << "Simulation took "<<delta.count()<< "s\n";
		std :: cout<<"Done!"<<std ::endl;
    		}

	MPI_Finalize();

	return 0;
}



