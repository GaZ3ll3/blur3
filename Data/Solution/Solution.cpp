/*
 * Solution.cpp
 *
 *  Created on: May 28, 2014
 *      Author: lurker
 */

#include "Solution.h"

Solution::Solution(SpatialMesh& Space, AngularMesh& Angle) {
	// TODO Auto-generated constructor stub
	flux.resize(Angle.ns);
	residue.resize(Angle.ns);
	for (register int i = 0; i < Angle.ns; i++){
		flux[i].resize(Space.nt);
		residue[i].resize(Space.nt);
		for (register int j = 0 ; j < Space.nt; j++){
			flux[i][j].resize(3);
			residue[i][j].resize(3);
		}
	}

}

Solution::~Solution() {
	// TODO Auto-generated destructor stub
}

