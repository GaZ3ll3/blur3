/*
 * SpatialMesh.cpp
 *
 *  Created on: May 26, 2014
 *      Author: lurker
 */

#include "SpatialMesh.h"

SpatialMesh::SpatialMesh(AngularMesh& Angle) {
	// TODO Auto-generated constructor stub
	std::fstream smeshfile("./smesh.txt", std::ios_base::in);
	double temp;
	smeshfile >> temp; nt = (int)temp;
	smeshfile >> temp; np = (int)temp;
	smeshfile >> temp; ne = (int)temp;
	so.resize(Angle.ns);
	for (int j = 0; j < Angle.ns; j++){
		so[j].resize(nt);
		for (int k = 0 ; k < nt; k++){
			smeshfile >> temp;
			so[j][k] = (int)temp - 1;
		} // end of loop k
	}// end of loop j
	p.resize(np);
	for (int j = 0; j < np; j++){
		p[j].resize(2);
		smeshfile >> p[j][0];
		smeshfile >> p[j][1];
	}// end of loop j
	t.resize(nt);
	for (int j = 0 ; j < nt; j++){
		t[j].resize(3);
		for (int k = 0; k < 3; k++){
			smeshfile >> temp;
			t[j][k] = (int)temp - 1;
		}// end of loop k
	}// end of loop j
	e.resize(ne);
	for (int j = 0; j < ne; j++){
		e[j].resize(4);
		for (int k = 1; k < 3; k++){
			smeshfile >> temp;
			e[j][k] = (int)temp - 1;
		}// end of loop k
	}//end of loop j
	smeshfile.close();
}

SpatialMesh::~SpatialMesh() {
	std::cout << "Destroyed Everything Spatial" << std::endl;
}

void SpatialMesh::Update(){
	c.resize(nt);
	for (int j = 0; j < nt; j++){
		c[j].resize(2);
		for (int k = 0; k < 2; k++){
			c[j][k] = (
					p[t[j][0]][k] +
					p[t[j][1]][k] +
					p[t[j][2]][k])/3.0;
		}
	}

	ec.resize(ne);
	for (int j = 0; j < ne; j++){
		ec[j].resize(2);
		for (int k = 0; k < 2; k++){
			ec[j][k] = (
					p[e[j][1]][k] +
					p[e[j][2]][k])/2.0;
		}
	}

	a.resize(nt);
	for (int j = 0; j < nt; j++){
		double x1 = p[t[j][0]][0]; double y1 = p[t[j][0]][1];
		double x2 = p[t[j][1]][0]; double y2 = p[t[j][1]][1];
		double x3 = p[t[j][2]][0]; double y3 = p[t[j][2]][1];
		a[j] = Area(x1,y1,x2,y2,x3,y3);
	}
}

void SpatialMesh::Update_adj(){
	IMATRIX temp_p2;
	temp_p2.resize(np);
	p2.resize(np);
	for (register int j = 0; j < np; j++){
		temp_p2[j].resize(TEMPSIZE);
	}

	for (register int j = 0; j < nt; j++){
		for (register int k = 0; k < 3; k++){
			temp_p2[t[j][k]][0] += 1; // each triangle shared by three nodes
			temp_p2[t[j][k]][temp_p2[t[j][k]][0]] = j;
		}
	}//end of loop k

	for (register int j = 0; j < np; j++){
		p2[j].resize(temp_p2[j][0]+1);
	}

	for (register int j =0; j < np; j++){
		for (register int k = 0; k < temp_p2[j][0]+1; k++){
			p2[j][k] = temp_p2[j][k];
		}
	}
}

void SpatialMesh::Update_boundary(){
//initialization
	e2.resize(ne);
	for (register int j = 0; j < ne; j++){
		e2[j].resize(2);
	}
	so2.resize(nt);
	for (register int j = 0; j < nt; j++){
		so2[j].resize(3);
		// initialized
		so2[j][0] = -1;
		so2[j][1] = -1;
		so2[j][2] = -1;
	}
	n.resize(ne);
	for (register int j = 0; j < ne; j++){
		n[j].resize(2);
	}
	ori.resize(ne);

	for (register int j = 0; j < ne; j++){
		int tri = locate_tri(p2[e[j][1]],p2[e[j][2]]);
		if (tri == -1){
			std::cout << "Not a boundary edge." << std::endl;
			break;
		}
        e[j][0]=tri;
        // local
        for (int k = 1; k < 3; k++){
        	for (int l = 0; l < 3; l++){
        		if (e[j][k] == t[tri][l]){
        			e2[j][k] = l;
        			break;
        		}
        	}
        }
        // inner node
        for (int l = 0; l < 3; l++){
        	if (t[tri][l]!= e[j][1] && t[tri][l]!=e[j][2]){
        		e[j][3]=t[tri][l];
        		so2[tri][l] = j;
        		break;
        	}
        }

        double nx = -(p[e[j][1]][1] - p[e[j][2]][1]);
        double ny = (p[e[j][1]][0] - p[e[j][2]][0]);

        ori[j] = 1;//counterclockwise

        if (nx * (p[e[j][3]][0] - p[e[j][2]][0]) + ny * (p[e[j][3]][1] - p[e[j][2]][1]) > 0){
        	nx = -nx;
        	ny = -ny;
        	ori = 0; // clockwise
        }
        n[j][0] = nx/sqrt(pow(nx,2) + pow(ny,2));
        n[j][1] = ny/sqrt(pow(nx,2) + pow(ny,2));
	}
}

void SpatialMesh::Update_boundary(AngularMesh& Angle){
	bd.resize(Angle.ns);
	for (register int j = 0; j < Angle.ns; j++){
		bd[j].resize(nt);
		for (register int k = 0; k < nt; k++){
			bd[j][k].resize(9);
			for (int l = 0 ; l < 9 ; l++){
				bd[j][k][l] = -1;
			}
		}
	}
	bd2.resize(Angle.ns);
	for (register int j = 0; j < Angle.ns; j++){
		bd2[j].resize(nt);
		for (register int k = 0; k < nt; k++){
			bd2[j][k].resize(3);
		}
	}
	for (register int j = 0 ; j < Angle.ns; j++){
		Flux(nt, Angle.a[j], p, p2, t, bd[j], bd2[j], so2);
	}
}

void SpatialMesh::Flux(int nt,DVECTOR& theta,DMATRIX& p,IMATRIX& p2,IMATRIX& t,IMATRIX& bd,DMATRIX& bd2,IMATRIX& so2){
	double a = theta[0];
	double b = theta[1];
	double dx, dy;

	for (register int j = 0; j < nt; j++){
		double x1 = p[t[j][0]][0]; double y1 = p[t[j][0]][1];
		double x2 = p[t[j][1]][0]; double y2 = p[t[j][1]][1];
		double x3 = p[t[j][2]][0]; double y3 = p[t[j][2]][1];

		// 1st edge "23"
		dx = y3 - y2; dy = x2 - x3;
		if(dx*(x1-x3)+dy*(y1-y3)>0){
			dx=-dx;
			dy=-dy;
		}
		bd2[j][0] = a*dx + b*dy;
		if (bd2[j][0] < 0 && so2[j][0] == -1){
			bd[j][0] = locate_tri(j, p2[t[j][1]], p2[t[j][2]]);
			if (bd[j][0]!=-1){
				//exists
				for(int k = 0; k < 3; k++){
					if(t[bd[j][0]][k]==t[j][1]){
						bd[j][1] = k;
						break;
					}
				}
				for(int k = 0; k < 3; k++){
					if(t[bd[j][0]][k]==t[j][2]){
						bd[j][2] = k;
						break;
					}
				}
			}
		}

		// 2nd edge "31"
		dx = y1 - y3; dy = x3 - x1;
		if(dx*(x2-x3)+dy*(y2-y3)>0){
			dx=-dx;
			dy=-dy;
		}
		bd2[j][1] = a*dx + b*dy;
		if (bd2[j][1] < 0 && so2[j][1] == -1){
			bd[j][3] = locate_tri(j, p2[t[j][0]], p2[t[j][2]]);
			if (bd[j][3]!=-1){
				//exists
				for(int k = 0; k < 3; k++){
					if(t[bd[j][3]][k]==t[j][2]){
						bd[j][4] = k;
						break;
					}
				}
				for(int k = 0; k < 3; k++){
					if(t[bd[j][3]][k]==t[j][0]){
						bd[j][5] = k;
						break;
					}
				}
			}
		}
		// 3rd edge "12"
		dx = y2 - y1; dy = x1 - x2;
		if(dx*(x3-x1)+dy*(y3-y1)>0){
			dx=-dx;
			dy=-dy;
		}
		bd2[j][2] = a*dx + b*dy;
		if (bd2[j][2] < 0 && so2[j][2] == -1){
			bd[j][6] = locate_tri(j, p2[t[j][1]], p2[t[j][2]]);
			if (bd[j][6]!=-1){
				//exists
				for(int k = 0; k < 3; k++){
					if(t[bd[j][6]][k]==t[j][0]){
						bd[j][7] = k;
						break;
					}
				}
				for(int k = 0; k < 3; k++){
					if(t[bd[j][6]][k]==t[j][1]){
						bd[j][8] = k;
						break;
					}
				}
			}
		}
	}
}

int SpatialMesh::locate_tri(int tri, IVECTOR& pt1, IVECTOR pt2){
	for (int i = 1; i < pt1[0] + 1; i++){
		for (int j = 1; j < pt2[0] + 1; j++){
			if (pt1[i] == pt2[j] && tri!=pt1[i]){
				return pt1[i];
			}
		}
	}
	return -1;
}

int SpatialMesh::locate_tri(IVECTOR& pt1, IVECTOR& pt2){
	for (int i = 1; i < pt1[0] + 1; i++){
		for (int j = 1; j < pt2[0] + 1; j++){
			if (pt1[i] == pt2[j]){
				return pt1[i];
			}
		}
	}
	return -1;
}

void SpatialMesh::Show(){

	std::cout << "np = " << np << std::endl;
	std::cout << "nt = " << nt << std::endl;
	std::cout << "ne = " << ne << std::endl;

	for (std::size_t j = 0 ; j < so.size(); j++){
		std::cout << "Sweeping order for angle[" << j << "]:" ;
		for (std::size_t k = 0; k < so[j].size(); k++){
			std::cout << so[j][k] << " ";
		}
		std::cout << std::endl;
	}

	for (std::size_t j = 0; j < c.size(); j++){
		std::cout << "centers: [" << j << "] = " << "(" << c[j][0] << "," << c[j][1] << ")"<< std::endl;
	}
	for (std::size_t j = 0; j < a.size(); j++){
		std::cout << "Area["<< j << "] = " << a[j] << std::endl;
	}
	for (std::size_t j = 0; j < p2.size(); j++){
		std::cout << "number: " << p2[j][0] << std::endl;
		for (std::size_t k = 1; k < p2[j].size(); k++){
			std::cout << "Triangles: " << p2[j][k] << std::endl;
		}
	}
	for (std::size_t j = 0; j < ec.size(); j++){
		for (std::size_t k = 0; k < ec[j].size(); k++){
			std::cout << "Edge centers: " << ec[j][k] << std::endl;
		}
	}
	for (std::size_t j = 0; j < ori.size(); j++){
		std::cout << "ori["<< j << "] = " << ori[j] << std::endl;
	}
	for (std::size_t j = 0; j < n.size(); j++){
		std::cout << "normals[" << j << "] = " << "(" << n[j][0] << "," << n[j][1] << ")"<< std::endl;
	}
	for (std::size_t j = 0; j < e.size(); j++){
		std::cout << "edges: [" << j << "] = " << "(" << e[j][0] << "," << e[j][1] << "," << e[j][2] << "," << e[j][3] << ")"<< std::endl;
	}
	for (std::size_t j = 0; j < so2.size(); j++){
		std::cout << "so2: [" << j << "] = " << so2[j][0] << "," << so2[j][1] << "," << so2[j][2] << std::endl;
	}
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

