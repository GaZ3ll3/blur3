/*
 * Solution.h
 *
 *  Created on: May 28, 2014
 *      Author: lurker
 */

#ifndef SOLUTION_H_
#define SOLUTION_H_

#include "../../Util.h"
#include "../../Mesh/Angular/src/AngularMesh.h"
#include "../../Mesh/Spatial/src/SpatialMesh.h"

class Solution {
public:
	Solution(SpatialMesh&, AngularMesh&);
	virtual ~Solution();
	DTENSOR sol;
	DTENSOR flux;
	DTENSOR residue;

};

#endif /* SOLUTION_H_ */
