#include "Common/Options.h"

namespace opt {
	/** Verbose output */
	int verbose;

	unsigned streakThreshold = 3;
	unsigned threads = 1;
	FilterType filterType = BLOOMFILTER;
	const ID EMPTY = 0;
	const ID COLLI = std::numeric_limits<ID>::max();
	std::vector<std::string> sseeds;
	unsigned allowMisses = 0;
	bool idByFile = false;
	bool colliIDs = false;
	bool colliAnalysis = false;
}
