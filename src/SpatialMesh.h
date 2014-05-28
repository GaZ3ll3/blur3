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
typedef struct spatialmesh{
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
    IMATRIX e2;
    IMATRIX so2;
    DMATRIX n;
    ITENSOR bd;
    ITENSOR bd2;
} spatialmesh;

class SpatialMesh {
public:
	SpatialMesh(int,int,AngularMesh);
	virtual ~SpatialMesh();
	void Show();
	void Update();
	void Mapping();
	void Mapping(spatialmesh&, spatialmesh&, IMATRIX&);
	void Mapping(spatialmesh&, spatialmesh&, IMATRIX&, DTUPLE&, DTUPLE&);
	void Edge();
	void Edge(int, int, IMATRIX&,IMATRIX&,DMATRIX&, IMATRIX&, IMATRIX&, IMATRIX&, DMATRIX&, IVECTOR&);
	std::size_t locate_min(DVECTOR&);
	double find_max(DVECTOR&);
	double Area(double x1, double y1, double x2, double y2, double x3, double y3);
	int slevel, slevel0, alevel, alevel0;
	int ds;
	std::vector<spatialmesh> smesh;
};



#endif /* SPATIALMESH_H_ */
