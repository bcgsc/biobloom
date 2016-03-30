#ifndef COMMON_OPTIONS_H
#define COMMON_OPTIONS_H 1

#include <stdint.h>
#include <limits>

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
}
#endif
