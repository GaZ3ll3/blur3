/*
 * main.cpp
 *
 *  Created on: May 28, 2014
 *      Author: lurker
 */
#include "Mesh/Angular/src/AngularMesh.h"
#include "Mesh/Spatial/src/SpatialMesh.h"

int main(){

	AngularMesh Angle;
	Angle.Show();
	SpatialMesh Space(Angle);
	Space.Update();
	Space.Update_adj();
	Space.Update_boundary();
	Space.Update_boundary(Angle);
	Space.Show();
	return 0;
}



