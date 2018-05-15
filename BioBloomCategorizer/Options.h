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
enum OutType {NONE ,FASTQ, FASTA, TSV};
extern bool inclusive;
extern double score;
extern std::string outputPrefix;
extern std::string filePostfix;
extern OutType outputType;
extern bool minHitOnly;
extern unsigned maxGroupSize;
extern int debug;
extern unsigned multiThresh;
extern bool inverse;

extern std::string filtersFile;
extern bool paired;
extern bool stdout;
}
#endif
