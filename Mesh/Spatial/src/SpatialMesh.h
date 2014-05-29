/*
 * SpatialMesh.h
 *
 *  Created on: May 26, 2014
 *      Author: lurker
 */

#ifndef SPATIALMESH_H_
#define SPATIALMESH_H_
#include "../../Angular/src/AngularMesh.h"
#include "../../../Util.h"

class SpatialMesh {
public:
	SpatialMesh(AngularMesh&);
	virtual ~SpatialMesh();
	void Show();
	void Update();
	void Update_adj();
	void Update_boundary();
	void Update_boundary(AngularMesh&);
	void Flux(int nt,DVECTOR& theta,DMATRIX& p,IMATRIX& p2,IMATRIX& t,IMATRIX& bd,DMATRIX& bd2,IMATRIX& so2);
	std::size_t locate_min(DVECTOR&);
	int locate_tri(IVECTOR& , IVECTOR&);
	int locate_tri(int tri, IVECTOR& pt1, IVECTOR pt2);
	double find_max(DVECTOR&);
	double Area(double x1, double y1, double x2, double y2, double x3, double y3);
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
//    IMATRIX smap;
//    DTUPLE cf;
//    DTUPLE fc;
    IVECTOR ori;
    IMATRIX e2;
    IMATRIX so2;
    DMATRIX n;
    ITENSOR bd;
    DTENSOR bd2;
};



#endif /* SPATIALMESH_H_ */
