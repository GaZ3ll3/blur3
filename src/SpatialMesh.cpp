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

void SpatialMesh::Show(){
	for (int i = 0; i <= slevel; i++)
	{
		std::cout << smesh[i].np << std::endl;
		std::cout << smesh[i].nt << std::endl;
		std::cout << smesh[i].ne << std::endl;
	}
}

void SpatialMesh::Update(){
/*
 *  ec[ne][2], a[nt], p2[np][p2[np][0] + 1]
 */
	int temp;
	for (int i = 0 ; i <= slevel; i++){
		temp = smesh[i].nt;
		smesh[i].c.resize(temp);
		for (int j = 0; j < temp; j++){
			smesh[i].c[j].resize(2);
			for (int k = 0; k < 2; k++){
				smesh[i].c[j][k] =(
						smesh[i].p[smesh[i].t[j][0]][k] +
						smesh[i].p[smesh[i].t[j][1]][k] +
						smesh[i].p[smesh[i].t[j][2]][k])/3.0;
			}
		}// end of loop j
		temp = smesh[i].ne;
		smesh[i].ec.resize(temp);
		for (int j  =0; j < temp; j++){
			smesh[i].ec[j].resize(2);
			for (int k = 0 ; k < 2; k++){
				smesh[i].ec[j][k] = (
						smesh[i].p[smesh[i].e[j][0]][k] +
						smesh[i].p[smesh[i].e[j][1]][k])/2.0;
			}
		}//end of loop j
		temp = smesh[i].nt;
		smesh[i].a.resize(temp);
		for (int j =0; j < temp; j++){
			double x1 = smesh[i].p[smesh[i].t[j][0]][0];
			double x2 = smesh[i].p[smesh[i].t[j][1]][0];
			double x3 = smesh[i].p[smesh[i].t[j][2]][0];
			double y1 = smesh[i].p[smesh[i].t[j][0]][1];
			double y2 = smesh[i].p[smesh[i].t[j][1]][1];
			double y3 = smesh[i].p[smesh[i].t[j][2]][1];
			smesh[i].a[j] = Area(x1,y1,x2,y2,x3,y3);
		}// end of loop j
		temp = smesh[i].np;
		smesh[i].p2.resize(temp);
		for (int j = 0; j < temp; j++){
			// dynamically allocate space, if failed, then spatial mesh is not acceptable. Angle too far away from Pi/3.
			smesh[i].p2[j].resize(TEMPSIZE);
		}
	}//end of loop i
}

double SpatialMesh::Area(double x1, double y1, double x2, double y2, double x3, double y3){
	return fabs(0.50*(x2*y3 + x3*y1 + x1*y2 - x2*y1 -x3*y2 - y3*x1));
}

