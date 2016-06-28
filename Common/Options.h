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
	extern int verbose;
	extern unsigned streakThreshold;
	extern unsigned threads;
	extern FilterType filterType;
	extern const ID EMPTY;
	extern const ID COLLI;
	extern std::vector<std::string> sseeds;
	extern unsigned allowMisses;

	//options of BBM only
	//TODO: move to own header file
	extern bool idByFile;
	extern bool colliIDs;
	extern bool colliAnalysis;
}
#endif
