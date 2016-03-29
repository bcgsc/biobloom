/*
 * Options.h
 *
 *  Created on: Mar 24, 2016
 *      Author: cjustin
 */

#ifndef OPTIONS_H
#define OPTIONS_H 1

#include <stdint.h>
#include <string>

/**
 * Global variables that are mostly constant for the duration of the
 * execution of the program.
 */
namespace opt {
extern double score;
extern std::string outputPrefix;
extern std::string filePostfix;
extern std::string outputType;
}
#endif
