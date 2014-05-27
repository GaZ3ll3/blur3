/*
 * SpatialMesh.cpp
 *
 *  Created on: May 26, 2014
 *      Author: lurker
 */

#include "SpatialMesh.h"

SpatialMesh::SpatialMesh(int sl, int sl0, AngularMesh Angle) {
	// TODO Auto-generated constructor stub
	slevel = sl; slevel0 = sl0;
	ds = slevel - slevel0;
	std::fstream smeshfile("src/smesh.txt", std::ios_base::in);
	for (int i = 0 ; i <= slevel; i++){
		struct spatialmesh refine; double temp; smeshfile >> temp;
		/* update np, ne, nt */
		refine.nt = (int)temp;smeshfile >> temp;
		refine.np = (int)temp;smeshfile >> temp;
		refine.ne = (int)temp;
		/* update so */
		refine.so.resize(Angle.amesh[Angle.alevel].ns);
		for (int j = 0; j < Angle.amesh[Angle.alevel].ns; j++){
			refine.so[j].resize(refine.nt);
			for (int k = 0 ; k < refine.nt; k++){
				smeshfile >> refine.so[j][k];
			} // end of loop k
		}// end of loop j
		/* update */





		smesh.push_back(refine);
	}// end of loop i

}

SpatialMesh::~SpatialMesh() {
	// TODO Auto-generated destructor stub
}

