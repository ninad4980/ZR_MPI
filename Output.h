
#ifndef _OUTPUT_H
#define _OUTPUT_H

#include "Domain.h"
#include "Field.h"
#include <cstring>
#include <string>
#include <vector>
#include <memory>

#include "Species.h"
using namespace std;

class Output 
{
public:
 	
	Output(Domain &domain) : domain(domain) {
	
	Nr = domain.Nr ;

	Nzee = domain.Nzee;
	Nz_phys = domain.Nz_phys;
	Nz_local = domain.Nz_local;
	}

	
	void write2file_Object();

	void write2file_SField2d(Field* F, char* fileName );	
	void write2file_SField3d(Field* F, char* fileName );

	void write2file_VField2d(FieldCyl* F, const char* fileName );

	/*saves simulation state data vs time*/
	void outputPart(
	    std::vector<std::unique_ptr<Species>> &speciesList,
	    const char* fileName
		);
	
	void outputPart2(
		KineticSpecies &species,
		const char* fileName
		);	
		
protected:
	Domain &domain;
	int Nr, Ntheta, Nzee, Nz_phys, Nz_local;
	
	std :: string path="results/";	
    
};

#endif
