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
	AngularMesh Angle(3,1);
//	Angle.show();
	SpatialMesh Space(0,0,Angle);
	Space.show();
	return 0;
}


