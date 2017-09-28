#ifndef COMMON_OPTIONS_H
#define COMMON_OPTIONS_H 1

#include <stdint.h>
#include <limits>
#include <string>
#include <vector>

enum FilterType {BLOOMFILTER, BLOOMMAP};

typedef uint16_t ID;

/**
 * Global variables that are mostly constant for the duration of the
 * execution of the program.
 */
namespace opt {
	extern bool colourSpace;
	extern int verbose;
	extern unsigned streakThreshold;
	extern unsigned threads;
	extern FilterType filterType;
	extern const ID EMPTY;
	extern const ID COLLI;
	extern std::vector<std::string> sseeds;
	extern unsigned allowMisses;
	extern double baitThreshold;
	extern unsigned progItrns;
	extern std::vector<std::string> fileList1;
	extern std::vector<std::string> fileList2;
	extern unsigned fileInterval;
	extern double fpr;
	extern bool noRep;
	
	//options of BBM only
	//TODO: move to own header file
	extern bool colliIDs;
	extern bool colliAnalysis;
	extern bool idByFile;


	//options for new refactored code
	extern std::string prefix;
	extern unsigned kmerSize;
	extern unsigned hashNum;
	extern unsigned numEle;
	extern std::string subtract;
	extern double pScore;
}

#endif
