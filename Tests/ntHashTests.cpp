/*
 *      Author: cjustin
 */

#include <string>
#include <assert.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include "btl_bloomfilter/MIBloomFilter.hpp"
#include "btl_bloomfilter/vendor/stHashIterator.hpp"
#include "Common/Options.h"
#if _OPENMP
# include <omp.h>
#endif

using namespace std;

int main()
{
	std::vector<string> ss;
	ss.push_back("111100001");
	ss.push_back("100001111");
	string seqFW = "TCAAATCTAA";
	string seqRV = "TTAGATTTGA";
	stHashIterator itr1(seqFW, MIBloomFilter<ID>::parseSeedString(ss), ss.size(), 1, ss[0].size());
	stHashIterator itr2(seqRV, MIBloomFilter<ID>::parseSeedString(ss), ss.size(), 1, ss[0].size());
	cerr << (*itr1)[0] << " " << (*itr2)[1] << endl;
	cerr << (*itr1)[1] << " " << (*itr2)[0] << endl;

	assert((*itr1)[0] != (*itr2)[1]);
	assert((*itr1)[1] != (*itr2)[0]);

	++itr1;
	cerr << (*itr1)[0] << " " << (*itr2)[1] << endl;
	cerr << (*itr1)[1] << " " << (*itr2)[0] << endl;

	assert((*itr1)[0] == (*itr2)[1]);
	assert((*itr1)[1] == (*itr2)[0]);
}
