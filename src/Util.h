/*
 * Util.h
 *
 *  Created on: May 27, 2014
 *      Author: lurker
 */

#ifndef UTIL_H_
#define UTIL_H_

#include <iostream>
#include <fstream>
#include <vector>
#include <valarray>

typedef std::valarray<std::valarray<int>> IMATRIX;
typedef std::valarray<std::valarray<double>> DMATRIX;
typedef std::valarray<double> DVECTOR;
typedef std::valarray<int> IVECTOR;
typedef std::valarray<std::valarray<std::valarray<std::valarray<double>>>> DTUPLE;
typedef std::valarray<std::valarray<std::valarray<int>>> ITENSOR;
typedef std::valarray<std::valarray<std::valarray<double>>> DTENSOR;

#define TEMPSIZE 20

#endif /* UTIL_H_ */
