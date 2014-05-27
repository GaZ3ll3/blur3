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
		struct spatialmesh refine;
		double temp; smeshfile >> temp;
		/* update np, ne, nt */
		refine.nt = (int)temp;smeshfile >> temp;
		refine.np = (int)temp;smeshfile >> temp;
		refine.ne = (int)temp;
		/* update so[ns][nt] */
		refine.so.resize(Angle.amesh[Angle.alevel].ns);
		for (int j = 0; j < Angle.amesh[Angle.alevel].ns; j++){
			refine.so[j].resize(refine.nt);
			for (int k = 0 ; k < refine.nt; k++){
				smeshfile >> temp;
				refine.so[j][k] = (int)temp;
			} // end of loop k
		}// end of loop j
		/* update p[np][2] */
		refine.p.resize(refine.np);
		for (int j = 0; j < refine.np; j++){
			refine.p[j].resize(2);
			smeshfile >> refine.p[j][0];
			smeshfile >> refine.p[j][1];
		}// end of loop j
		/*
		 * update t[nt][3]
		 */
		refine.t.resize(refine.nt);
		for (int j = 0 ; j < refine.nt; j++){
			refine.t[j].resize(3);
			for (int k = 0; k < 3; k++){
				smeshfile >> temp;
				refine.t[j][k] = (int)temp;
			}// end of loop k
		}// end of loop j
		/*
		 * update e[ne][4]
		 */
		refine.e.resize(refine.ne);
		for (int j = 0; j < refine.ne; j++){
			refine.e.resize(2);
			for (int k = 0; k < 2; k++){
				smeshfile >> temp;
				refine.t[j][k] = (int)temp;
			}// end of loop k
		}//end of loop j
		smesh.push_back(refine);
	}// end of loop i
	smeshfile.close();
}

SpatialMesh::~SpatialMesh() {
	// TODO Auto-generated destructor stub
}

void SpatialMesh::show(){
	for (int i = 0; i <= slevel; i++)
	{
		std::cout << smesh[i].np << std::endl;
		std::cout << smesh[i].nt << std::endl;
		std::cout << smesh[i].ne << std::endl;
	}
}

