/*
 * Source.cpp
 *
 *  Created on: May 28, 2014
 *      Author: lurker
 */

#include "Source.h"

// Right now only supports isotropic point sources

Source::Source(SpatialMesh& Space, AngularMesh& Angle) {
	// TODO Auto-generated constructor stub
	RHS.resize(Angle.ns);
	q.resize(Angle.ns);
	for (register int i = 0; i < Angle.ns; i++){
		RHS[i].resize(Space.nt);
		q[i].resize(Space.ne);
		for (register int j = 0; j < Space.nt; j++){
			RHS[i][j].resize(3);
		}
		for (register int j = 0; j < Space.ne; j++){
			q[i][j].resize(2);
		}
	}

	// Assembly light source
	std::fstream sourcefile("./source.txt", std::ios_base::in);
	sourcefile >> lumin.type;
}

Source::~Source() {
	// TODO Auto-generated destructor stub
}

double Source::distance(double& x1, double& y1, double& x2, double& y2){
	return sqrt(  (x1 - x2)*(x1 - x2) + ( y1 - y2 ) * (y1 - y2) );
}
