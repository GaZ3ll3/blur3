/*
 * AngularMesh.h
 *
 *  Created on: May 26, 2014
 *      Author: lurker
 */


#ifndef ANGULARMESH_H_
#define ANGULARMESH_H_

#include <iostream>
#include <fstream>
#include <vector>
#include <valarray>



struct angularmesh
{
	int ns;
	std::valarray<std::valarray<double > > a;
	std::valarray<std::valarray<double > > w;
};


class AngularMesh {
public:
	AngularMesh(int, int);
	virtual ~AngularMesh();
	std::vector<struct angularmesh> amesh;
	int da;
	int alevel;
	int alevel0;
	void show();

};

#endif /* ANGULARMESH_H_ */
