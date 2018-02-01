/*
 * Options.cpp
 *
 *  Created on: Mar 24, 2016
 *      Author: cjustin
 */

#include <Options.h>

using namespace std;
#include <limits>

namespace opt {
bool inclusive = false;
double score = 0.15;
string outputPrefix = "";
string filePostfix = "";
string outputType = "";
bool minHitOnly = false;
unsigned minHit = 0;
unsigned delta = 0;
unsigned maxGroupSize = std::numeric_limits<unsigned>::max();
bool inverse = false;
}


