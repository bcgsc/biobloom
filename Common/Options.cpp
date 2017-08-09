#include "Common/Options.h"

namespace opt {
	/** Colour space sequences */
	bool colourSpace;

	/** Verbose output */
	int verbose;

	unsigned streakThreshold = 3;

	unsigned threads = 1;

	double baitThreshold = -1;
	unsigned progItrns = 1;

	std::vector<std::string>fileList1;
	std::vector<std::string>fileList2;

	unsigned fileInterval = 10000000;
	double fpr = 0.0075;
	bool noRep = false;
}
