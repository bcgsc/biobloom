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

/** for modes of filtering */
enum FilteringMode {
	STD, ORDERED, BESTHIT, SCORES
};

extern bool inclusive;
extern double score;
extern std::string outputPrefix;
extern std::string filePostfix;
extern OutType outputType;
extern FilteringMode mode;
extern bool minHitOnly;
extern unsigned maxGroupSize;
extern int debug;
extern double multiThresh;
extern bool inverse;

extern double minFPR;

extern std::string filtersFile;
extern bool paired;
extern bool stdout;
extern bool bestHitCountAgree;

extern unsigned minCountNonSatCount;
extern unsigned frameMatches;

extern bool hitOnly;
}
#endif
