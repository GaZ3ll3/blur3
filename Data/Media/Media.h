/*
 * Media.h
 *
 *  Created on: May 28, 2014
 *      Author: lurker
 */

#ifndef MEDIA_H_
#define MEDIA_H_

#include "../../Util.h"
#include "../../Mesh/Angular/src/AngularMesh.h"
#include "../../Mesh/Spatial/src/SpatialMesh.h"

class Media {
public:
	Media(SpatialMesh&, AngularMesh&);
	virtual ~Media();
	DMATRIX ua;
	DMATRIX us;
	void Show();

};

#endif /* MEDIA_H_ */
