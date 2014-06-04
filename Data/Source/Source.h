/*
 * Source.h
 *
 *  Created on: May 28, 2014
 *      Author: lurker
 */

#ifndef SOURCE_H_
#define SOURCE_H_

#include "../../Util.h"
#include "../../Mesh/Angular/src/AngularMesh.h"
#include "../../Mesh/Spatial/src/SpatialMesh.h"

typedef struct source{
		int type;
		int n;
		DVECTOR a;
		DVECTOR t;
		DMATRIX i;
} source;


class Source {
public:
	Source(SpatialMesh&, AngularMesh&);
	virtual ~Source();
	DTENSOR RHS;
	DTENSOR q;
	source  lumin;
	source  boundary;
	double distance(double&,double&,double&,double&);

};

#endif /* SOURCE_H_ */
