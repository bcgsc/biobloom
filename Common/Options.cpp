#include "Common/Options.h"

namespace opt {
	/** Colour space sequences */
	bool colourSpace;

	/** Verbose output */
	int verbose;

	unsigned streakThreshold = 3;
	unsigned threads = 1;

	FilterType filterType = BLOOMFILTER;
	const ID EMPTY = 0;
	const ID COLLI = std::numeric_limits<ID>::max();
	std::vector<std::string> sseeds;
	unsigned allowMisses = 0;

	IDType idType = HEADER;
	bool colliIDs = false;
	bool colliAnalysis = false;
	bool idByFile = false;

	double baitThreshold = -1;
	unsigned progItrns = 1;

	std::vector<std::string>fileList1;
	std::vector<std::string>fileList2;

	unsigned fileInterval = 10000000;
	double fpr = 0.0078125;
	bool noRep = false;
	
	std::string prefix = "";
	unsigned kmerSize = 25;
	unsigned hashNum = 0;
	unsigned numEle = 0;
	std::string subtract = "";
}
