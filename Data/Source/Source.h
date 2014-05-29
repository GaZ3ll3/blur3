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
class Source {
public:
	Source(SpatialMesh&, AngularMesh&);
	virtual ~Source();
	DTENSOR RHS;
	DTENSOR q;
};

#endif /* SOURCE_H_ */
