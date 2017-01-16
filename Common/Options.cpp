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

	bool fastIO = true;

	double baitThreshold = -1;
	unsigned progItrns = 1;

	double fpr = 0.0078125;
	std::string prefix = "";
	unsigned kmerSize = 64;
	unsigned hashNum = 0;
	unsigned numEle = 0;
	std::string subtract = "";

	double pScore = 0.5;
	SeqEval::EvalMode mode = SeqEval::EVAL_STANDARD;
}
