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
		spatialmesh refine;
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
				refine.so[j][k] = (int)temp - 1;
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
				refine.t[j][k] = (int)temp - 1;
			}// end of loop k
		}// end of loop j
		/*
		 * update e[ne][4]
		 */
		refine.e.resize(refine.ne);
		for (int j = 0; j < refine.ne; j++){
			refine.e[j].resize(2);
			for (int k = 0; k < 2; k++){
				smeshfile >> temp;
				refine.e[j][k] = (int)temp - 1;
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
		std::cout << i<<"th level:" << std::endl;
		std::cout << "np = " << smesh[i].np << std::endl;
		std::cout << "nt = " << smesh[i].nt << std::endl;
		std::cout << "ne = " << smesh[i].ne << std::endl;

		for (std::size_t j = 0 ; j < smesh[i].so.size(); j++){
			std::cout << "Sweeping order for angle[" << j << "]:" ;
			for (std::size_t k = 0; k < smesh[i].so[j].size(); k++){
				std::cout << smesh[i].so[j][k] << " ";
			}
			std::cout << std::endl;
		}

		for (std::size_t j = 0; j < smesh[i].c.size(); j++){
			std::cout << "centers: [" << j << "] = " << "(" << smesh[i].c[j][0] << smesh[i].c[j][1] << ")"<< std::endl;
		}
		for (std::size_t j = 0; j < smesh[i].a.size(); j++){
			std::cout << "Area["<< j << "] = " << smesh[i].a[j] << std::endl;
		}
		for (std::size_t j = 0; j < smesh[i].p2.size(); j++){
			std::cout << "number: " << smesh[i].p2[j][0] << std::endl;
			for (std::size_t k = 1; k < smesh[i].p2[j].size(); k++){
				std::cout << "Triangles: " << smesh[i].p2[j][k] << std::endl;
			}
		}
		for (std::size_t j =0 ; j < smesh[i].smap.size(); j++){
			std::cout << "smap: " << smesh[i].smap[j][0] << std::endl;
			for (std::size_t k = 1; k < smesh[i].smap[j].size(); k++){
				std::cout << "smap[" << j << "][" << k << "] = " << smesh[i].smap[j][k] << std::endl;
			}
		}
	}
}

void SpatialMesh::Update(){
/*
 * update some information
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
		// Auxiliary matrix, will free at last
		IMATRIX p2;
		p2.resize(temp);
		smesh[i].p2.resize(temp);
		for (int j = 0; j < temp; j++){
			// dynamically allocate space, if failed, then spatial mesh is not acceptable. Angle too far away from Pi/3.
			p2[j].resize(TEMPSIZE);
		}
		for (int k = 0; k < smesh[i].nt; k++){
			for (int l = 0; l < 3; l++){
				// traversal of triangle
				p2[smesh[i].t[k][l]][0] += 1; // each triangle shared by three nodes
				p2[smesh[i].t[k][l]][p2[smesh[i].t[k][l]][0]] = k;
//				std::cout << "p2[" << smesh[i].t[k][l] << "][0] = " << p2[smesh[i].t[k][l]][0]  << std::endl;
//				std::cout << "p2[" << smesh[i].t[k][l] << "][" << p2[smesh[i].t[k][l]][0] << "] = " <<  p2[smesh[i].t[k][l]][p2[smesh[i].t[k][l]][0]] << std::endl;
			}
		}//end of loop k
		for (int j = 0; j < temp; j++){
			smesh[i].p2[j].resize(p2[j][0]+1);
		}
		for (int j =0; j < temp; j++){
			for (int k = 0; k < p2[j][0]+1; k++){
				smesh[i].p2[j][k] = p2[j][k];
			}
		}
	}//end of loop i
}

void SpatialMesh::Mapping(){
	// initialization
	IMATRIX smap;
	for (int i = 1; i <= slevel; i++){
		smap.resize(smesh[i-1].nt);
		smesh[i].smap.resize(smesh[i-1].nt);
		for (int j = 0; j < smesh[i-1].nt ; j++){
			smap[j].resize(TEMPSIZE);
		}// end of loop j

		Mapping(smesh[i-1], smesh[i], smap);

		for (int j = 0; j < smesh[i-1].nt; j++){
			smesh[i].smap[j].resize(smap[j][0]+1);
			for (int k = 0; k < smap[j][0] + 1; k++){
				smesh[i].smap[j][k] = smap[j][k];
			}// end of loop k
		}//end of loop j

		smesh[i].cf.resize(smesh[i-1].nt);

		for (int j = 0; j < smesh[i-1].nt; j++){
			smesh[i].cf[j].resize(3);
			for (int k =0 ; k < 3; k++){
				smesh[i].cf[j][k].resize(smesh[i].smap[j][0]);
				for (int l=0; l < smesh[i].smap[j][0]; l++){
					smesh[i].cf[j][k][l].resize(3);
				}// end of loop l
			}//end of loop k
		}//end of loop j
		smesh[i].fc.resize(smesh[i-1].nt);
		for (int j = 0; j < smesh[i-1].nt; j++){
			smesh[i].fc[j].resize(3);
			for (int k =0 ; k < 3; k++){
				smesh[i].fc[j][k].resize(smesh[i].smap[j][0]);
				for (int l=0; l < smesh[i].smap[j][0]; l++){
					smesh[i].fc[j][k][l].resize(3);
				}//end of loop l
			}//end of loop k
		}//end of loop j

		Mapping(smesh[i-1], smesh[i], smesh[i].smap, smesh[i].cf, smesh[i].fc);

	}// end of loop i
}

void SpatialMesh::Mapping(spatialmesh& cmesh, spatialmesh& fmesh, IMATRIX& smap){
	// initialization
	int nt_c = cmesh.nt;
	int tri;
	double dmax, x1,x2,x3,y1,y2,y3, area1,area2,area3;
	DVECTOR distance; distance.resize(nt_c);
	for (int i = 0; i < fmesh.nt; i++){
		int flag = 0;
		double x = fmesh.c[i][0];
		double y = fmesh.c[i][1];
		for (int j = 0; j < nt_c; j++){
			distance[j] = pow(x - cmesh.c[j][0],2) + pow(y - cmesh.c[j][1],2);
		}
		// linear time
		dmax = find_max(distance);
		// linear search maximum will be faster, since it is neighborhood related only. Finished in linear time.
		// quick sort needs NlogN time complexity.
		for (int j = 0; j < nt_c ; j++){
			tri = locate_min(distance);
			distance[(int)tri] = dmax;
			x1 = cmesh.p[cmesh.t[tri][0]][0]; y1 = cmesh.p[cmesh.t[tri][0]][1];
			x2 = cmesh.p[cmesh.t[tri][1]][0]; y2 = cmesh.p[cmesh.t[tri][1]][1];
			x3 = cmesh.p[cmesh.t[tri][2]][0]; y3 = cmesh.p[cmesh.t[tri][2]][1];
			area1 = Area(x,y,x1,y1,x2,y2);
			area2 = Area(x,y,x2,y2,x3,y3);
			area3 = Area(x,y,x3,y3,x1,y1);
			if ((cmesh.a[tri] - area1 - area2 -area3)/cmesh.a[tri] < TOL){
				flag = 1;
				smap[tri][0]+=1;
				smap[tri][smap[tri][0]]=i;
				break;
			}
		}
		if (flag == 0){
			std::cout << "mapping does not fit" << std::endl;
		}
	}
}

void SpatialMesh::Mapping(spatialmesh& cmesh, spatialmesh& fmesh, IMATRIX& smap, DTUPLE& cf, DTUPLE& fc){

}

void SpatialMesh::Edge(){

}
std::size_t SpatialMesh::locate_min(DVECTOR& vec){
	double dmin = vec[0]; std::size_t ind = 0;
	for (std::size_t i = 1; i < vec.size(); i++){
		if (dmin > vec[i]){
			dmin = vec[i];
			ind = i;
		}
	}
	return ind;
}

double SpatialMesh::find_max(DVECTOR& vec){
	double dmax = vec[0];
	for (std::size_t i = 1; i < vec.size(); i++){
		if (dmax < vec[i]){
			dmax = vec[i];
		}
	}
	return dmax;
}


double SpatialMesh::Area(double x1, double y1, double x2, double y2, double x3, double y3){
	return fabs(0.50*(x2*y3 + x3*y1 + x1*y2 - x2*y1 -x3*y2 - y3*x1));
}

