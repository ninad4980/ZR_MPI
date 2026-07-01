/*defines the simulation domain*/
#ifndef _GEOMETRY_H
#define _GEOMETRY_H
#include <string>
#include <vector>
#include <list>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <fstream>
#include "Domain.h"
#include "Field.h"
using namespace Const;		
/*node type*/
enum NodeType {NORMAL,OPEN};

/*definition of a node*/
struct Node
{
	Node(double x, double y, double z) {pos[0]=x;pos[1]=y;pos[2]=z;type=NORMAL;}	
	double pos[3];	/*node position*/
	NodeType type;
	double volume;	/*node volume*/
};

/*definition of a tetrahedron*/


class GeometryManual
{

public:
	GeometryManual(Domain &domain) : domain(domain) {
	
	Nr = domain.Nr ;
	Ntheta = domain.Ntheta ; 
	Nz_phys = domain.Nz_phys;

	}
	void plasmaChamber();

	void AddObject_Inlet();

	void AddObject_electrode(
    double posE,
    double dE,
    double rAperture,
    double rOut,
    double phi_electrode,
    ObjectType type
	);

protected:
	Domain &domain;
	int Nr, Ntheta, Nz_phys;	
};























#endif
