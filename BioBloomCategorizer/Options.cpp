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
OutType outputType = NONE;

FilteringMode mode = STD;

bool minHitOnly = false;
unsigned maxGroupSize = std::numeric_limits<unsigned>::max();
int debug = 0;
double multiThresh = 1;
bool inverse = false;

std::string filtersFile = "";
bool paired = false;
bool stdout = false;
bool bestHitCountAgree = false;

unsigned minCountNonSatCount = 0;
unsigned frameMatches = 1;

bool hitOnly = false;
}


