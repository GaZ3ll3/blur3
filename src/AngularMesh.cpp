/*
 * AngularMesh.cpp
 *
 *  Created on: May 26, 2014
 *      Author: lurker
 */

#include "AngularMesh.h"

AngularMesh::AngularMesh(int al, int al0) {
	// TODO Auto-generated constructor stub
	alevel = al;alevel0 = al0;
	da = alevel - alevel0;
	std::fstream ameshfile("src/amesh.txt", std::ios_base::in);
	// time complexity O(n)
	for (int i = 0 ; i <= alevel; i++){
		angularmesh refine;
		register double temp;
		ameshfile >> temp;
		refine.ns = (int)temp;
		refine.a.resize(refine.ns);
		for (int j = 0 ; j < refine.ns ; j++){
			refine.a[j].resize(3);
			for (int k = 0 ; k < 3; k++){
				ameshfile >> refine.a[j][k];
			}
		}
		refine.w.resize(refine.ns);
		for (int j = 0 ; j < refine.ns ; j++){
			refine.w[j].resize(refine.ns);
			for (int k = 0 ; k < refine.ns; k++){
				ameshfile >> refine.w[j][k];
			}
		}
		amesh.push_back(refine);
	}
	ameshfile.close();
}

AngularMesh::~AngularMesh() {
	std::cout << "Destroyed Everything Angular" << std::endl;
}


void AngularMesh::show()
{
	for (register int i = 0; i <= alevel; i++)
	{
		std::cout << i <<"th level grid of angles:" << amesh[i].ns << std::endl;
	}
}

