/*
 * SpatialMesh.h
 *
 *  Created on: May 26, 2014
 *      Author: lurker
 */

#ifndef SPATIALMESH_H_
#define SPATIALMESH_H_
#include "AngularMesh.h"
#include "Util.h"
struct spatialmesh{
    int nt;
    int np;
    int ne;
    IMATRIX so;
    DMATRIX p;
    IMATRIX t;
    IMATRIX e;
    DMATRIX c;
    DMATRIX ec;
    DVECTOR a;
    IMATRIX p2;
    IMATRIX smap;
    DTUPLE cf;
    DTUPLE fc;
    IVECTOR ori;
    ITENSOR e2;
    ITENSOR so2;
    DMATRIX n;
    ITENSOR bd;
    ITENSOR bd2;
};

class SpatialMesh {
public:
	SpatialMesh(int,int,AngularMesh);
	virtual ~SpatialMesh();
	int slevel, slevel0, alevel, alevel0;
	int ds;
	std::vector<struct spatialmesh> smesh;
};

#endif /* SPATIALMESH_H_ */
