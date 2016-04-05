#ifndef COMMON_OPTIONS_H
#define COMMON_OPTIONS_H 1

#include <stdint.h>
#include <limits>
#include <string>
#include <vector>

enum FilterType {BLOOMFILTER, BLOOMMAP};

typedef uint16_t ID;
typedef uint32_t PairID;
static const unsigned ID_BITS = 16;

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
}
#endif
