/*
 * main.cpp
 *
 *  Created on: May 26, 2014
 *      Author: lurker
 */

#include "Solver.h"
#include "AngularMesh.h"
#include "SpatialMesh.h"

int main(int argc, char** argv)
{
	std::cout << "Blur 2" << std::endl;
	/*
	 * Initialization of spatial mesh and angular mesh.
	 */
	AngularMesh Angle(3,1);
	Angle.show();
	SpatialMesh Space(3,3,Angle);
	Space.Update();
	Space.Mapping();
	Space.Edge();
	Space.Boundary(Angle);
	Space.Show();
	/*
	 * End of Initialization
	 */
	return 0;
}


