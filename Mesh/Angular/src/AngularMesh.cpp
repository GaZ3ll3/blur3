/*
 * AngularMesh.cpp
 *
 *  Created on: May 26, 2014
 *      Author: lurker
 */

#include "AngularMesh.h"

AngularMesh::AngularMesh() {
	// TODO Auto-generated constructor stub
	std::fstream ameshfile("./amesh.txt", std::ios_base::in);
	// time complexity O(n)

	register double temp;
	ameshfile >> temp;
	ns = (int)temp;
	a.resize(ns);
	for (int j = 0 ; j < ns ; j++){
		a[j].resize(3);
		for (int k = 0 ; k < 3; k++){
			ameshfile >> a[j][k];
		}
	}
	w.resize(ns);
	for (int j = 0 ; j < ns ; j++){
		w[j].resize(ns);
		for (int k = 0 ; k < ns; k++){
			ameshfile >> w[j][k];
		}
	}
	ameshfile.close();
}

AngularMesh::~AngularMesh() {
	std::cout << "Destroyed Everything Angular" << std::endl;
}


void AngularMesh::Show()
{
	std::cout << "grid of angles:" << ns << std::endl;
}

