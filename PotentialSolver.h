#ifndef _SOLVER_H
#define _SOLVER_H
#include <vector>
#include <memory>
#include "Domain.h"
#include "Species.h"
using namespace Const;
class PotentialSolver 
{
public:
	/*constructor, sets world*/
	PotentialSolver(Domain &domain, int max_it, double tol): 
		domain(domain), max_solver_it(max_it), tolerance(tol) {}
	
	
	void calcRhoC(
	    std::vector<std::unique_ptr<Species>> &speciesList);
	/*solves potential using Gauss-Seidel*/
	bool solvePhi();
	/*computes electric field = -gradient(phi)*/
	void computeEF();

protected:
	Domain &domain;
	unsigned max_solver_it;	//maximum number of solver iterations
	double tolerance;		//solver tolerance

	int Nr, Ntheta, Nzee;
	double r, theta;
	double Aijk, Arpjk,  Armjk, Aitk,  Aijz;
	double dr, dtheta, dzee;
	double phiNew;

	bool converged;
	double w ;
private:

};
#endif
