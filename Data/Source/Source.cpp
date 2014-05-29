/*
 * Source.cpp
 *
 *  Created on: May 28, 2014
 *      Author: lurker
 */

#include "Source.h"


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
}

Source::~Source() {
	// TODO Auto-generated destructor stub
}

