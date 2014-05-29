/*
 * AngularMesh.h
 *
 *  Created on: May 26, 2014
 *      Author: lurker
 */


#ifndef ANGULARMESH_H_
#define ANGULARMESH_H_

#include "../../../Util.h"

class AngularMesh {
public:
	AngularMesh();
	virtual ~AngularMesh();
	int ns;
	DMATRIX a;
	DMATRIX w;
	void Show();

};

#endif /* ANGULARMESH_H_ */
