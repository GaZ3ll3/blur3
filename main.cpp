/*
 * main.cpp
 *
 *  Created on: May 28, 2014
 *      Author: lurker
 */
#include "Mesh/Angular/src/AngularMesh.h"
#include "Mesh/Spatial/src/SpatialMesh.h"
#include "Data/Media/Media.h"
#include "Data/Source/Source.h"
#include "Data/Solution/Solution.h"
int main(){

	AngularMesh Angle;
	Angle.Show();
	SpatialMesh Space(Angle);
	Space.Update();
	Space.Update_adj();
	Space.Update_boundary();
	Space.Update_boundary(Angle);
	Space.Show();
	Media Median(Space, Angle);
	Median.Show();
	return 0;
}



