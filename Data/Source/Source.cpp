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
	// Variables
	double temp;
	DMATRIX location;
	DVECTOR intensity;

	// Main body
	std::fstream sourcefile("./source.txt", std::ios_base::in);
	sourcefile >> temp;
	lumin.type = (int) temp;
	if (lumin.type != 0){
		sourcefile >> temp;
		lumin.n = (int) temp;
		lumin.i.resize(lumin.n);
		if (lumin.type == 3){
			// isotropic nodal source everywhere
			for (int j = 0; j < lumin.n; j++ ){
				lumin.i[j].resize(3);
				for ( int k = 0 ; k < 3; k++){
					sourcefile >> lumin.i[j][k];
				}
			}
			sourcefile.close();
		}
		else { // type is 1 or 2
			if (lumin.type == 1){
				// point source with deterministic angle
				lumin.a.resize(lumin.n);
				for (int j = 0; j < lumin.n; j++){
					sourcefile >> temp;
					lumin.a[j] = (int) temp - 1;
				}
			}

			// default setting
			location.resize(lumin.n);
			for (int j = 0 ; j < lumin.n ; j++){
				location[j].resize(2);
				for (int k = 0 ; k < 2; k++){
					sourcefile >> location[j][k];
				}
			}

			intensity.resize(lumin.n);
			for (int j = 0 ; j < lumin.n; j++){
				sourcefile >> intensity[j];
			}

			// assign t and i


		}
	}


	std::fstream bsourcefile("./bsource.txt", std::ios_base::in);
	bsourcefile >> temp;
	bound.type = (int) temp;
	if (bound.type != 0){
		bsourcefile >> temp;
		bound.n = (int) temp;
		bound.i.resize(bound.n);
		for (int j = 0; j < bound.n; j++ ){
			bound.i[j].resize(2);
			for ( int k = 0 ; k < 2; k++){
				bsourcefile >> bound.i[j][k];
			}
		}
	}
	sourcefile.close();
}

Source::~Source() {
	// TODO Auto-generated destructor stub
}

double Source::distance(double& x1, double& y1, double& x2, double& y2){
	return sqrt(  (x1 - x2)*(x1 - x2) + ( y1 - y2 ) * (y1 - y2) );
}
