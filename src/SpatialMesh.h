/*
 * SpatialMesh.h
 *
 *  Created on: May 26, 2014
 *      Author: lurker
 */

#ifndef SPATIALMESH_H_
#define SPATIALMESH_H_

class SpatialMesh {
public:
	SpatialMesh(int,int);
	virtual ~SpatialMesh();
	int slevel, slevel0;
};

#endif /* SPATIALMESH_H_ */
