#include <vector>
#include "btl_bloomfilter/stHashIterator.hpp"
#include <string>
#include <iostream>
#include <cassert>

using namespace std;

int main()
{
	std::vector<string> ss;
	ss.push_back("111100001");
	ss.push_back("100001111");
	string seqFW = "TCAAATCTAA";
	string seqRV = "TTAGATTTGA";

	stHashIterator itr1(seqFW, parseSeed(ss), ss.size(), ss[0].size());
	stHashIterator itr2(seqRV, parseSeed(ss), ss.size(), ss[0].size());
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
