#ifndef COMMON_OPTIONS_H
#define COMMON_OPTIONS_H 1

#include <vector>
#include <string>

/**
 * Global variables that are mostly constant for the duration of the
 * execution of the program.
 */
namespace opt {
	extern bool colourSpace;
	extern int verbose;
	extern unsigned streakThreshold;
	extern unsigned threads;
	extern int rank;
	extern double baitThreshold;
	extern unsigned progItrns;
	extern std::vector<std::string> fileList1;
	extern std::vector<std::string> fileList2;
	extern unsigned fileInterval;
	extern double fpr;
	extern bool noRep;
}

#endif
