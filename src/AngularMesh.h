/*
 * AngularMesh.h
 *
 *  Created on: May 26, 2014
 *      Author: lurker
 */


#ifndef ANGULARMESH_H_
#define ANGULARMESH_H_

#include "Util.h"

typedef struct angularmesh
{
	int ns;
	DMATRIX a;
	DMATRIX w;
} angularmesh;


class AngularMesh {
public:
	AngularMesh(int, int);
	virtual ~AngularMesh();
	std::vector<angularmesh> amesh;
	int da;
	int alevel;
	int alevel0;
	void show();

};

#endif /* ANGULARMESH_H_ */
