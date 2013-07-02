/*
 * BloomFilterTests.cpp
 * Unit Tests for hashmanager and bloomfilter classes
 *  Created on: Aug 14, 2012
 *      Author: cjustin
 */

//#include "BloomFilter.h"
#include "Common/HashManager.h"
#include "Common/BloomFilter.h"
#include <string>
#include <assert.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

//returns memory of program in kb
int memory_usage() {
	int mem = 0;
	ifstream proc("/proc/self/status");
	string s;
	while (getline(proc, s), !proc.fail()) {
		if (s.substr(0, 6) == "VmSize") {
			stringstream convert(
					s.substr(s.find_last_of('\t'), s.find_last_of('k') - 1));
			if (!(convert >> mem)) {
				return 0;
			}
			return mem;
		}
	}
	return mem;
}

int main(int argc, char **argv) {
	HashManager hashMan;
	for (unsigned int i = 0; i < 12; ++i) {
		hashMan.addHashFunction("CityHash64", i);
	}
	const vector<size_t> &hashedValues = hashMan.multiHash("ATCGGGTCATCAACCAATAT");
	cout << hashedValues.size() << endl;
	for (unsigned int i = 0; i < hashedValues.size() - 1; ++i) {
		for (unsigned int j = i + 1; j < hashedValues.size(); ++j) {
			assert(hashedValues.at(j) != hashedValues.at(i));
		}
	}

	// end of testing hash manager, numbers should all be different.
	cout << "hashing tests done" << endl;

	//memory usage from before
	int memUsage = memory_usage();

	size_t filterSize = 1000000000;
	BloomFilter filter(filterSize, hashMan);
	filter.insert("ATCGGGTCATCAACCAATAT");
	filter.insert("ATCGGGTCATCAACCAATAC");
	filter.insert("ATCGGGTCATCAACCAATAG");
	filter.insert("ATCGGGTCATCAACCAATAA");

	//Check if filter is able to report expected results
	assert(filter.contains("ATCGGGTCATCAACCAATAT"));
	assert(filter.contains("ATCGGGTCATCAACCAATAC"));
	assert(filter.contains("ATCGGGTCATCAACCAATAG"));
	assert(filter.contains("ATCGGGTCATCAACCAATAA"));

	assert(!filter.contains("ATCGGGTCATCAACCAATTA"));
	assert(!filter.contains("ATCGGGTCATCAACCAATTC"));

	//should be size of bf (amortized)
	cout << memory_usage() - memUsage << "kb" << endl;

	cout << "de novo bf tests done" << endl;

	//Check storage can occur properly
	string filename = "/home/cjustin/workspace/TestData/bloomFilter.bf";
	filter.storeFilter(filename);
	ifstream ifile(filename.c_str());
	assert(ifile.is_open());
	ifile.seekg(0, ios::end); // move to end of file
	size_t fileSize = ifile.tellg(); // file size in bytes
	//file size should be same as filter size (Round to block size)
	if (filterSize % 64 > 0) {
		assert((filterSize + (64 - (filterSize% 64))) == fileSize*8);
	} else {
		assert(filterSize == fileSize*8);
	}
	ifile.close();

	//should be roughly same size still (amortized)
	cout << memory_usage() - memUsage << "kb" << endl;

	//check loading of stored filter
	BloomFilter filter2(filterSize, filename, hashMan);

	//should be double size of bf (amortized)
	cout << memory_usage() - memUsage << "kb" << endl;

	//Check if loaded filter is able to report expected results
	assert(filter2.contains("ATCGGGTCATCAACCAATAT"));
	assert(filter2.contains("ATCGGGTCATCAACCAATAC"));
	assert(filter2.contains("ATCGGGTCATCAACCAATAG"));
	assert(filter2.contains("ATCGGGTCATCAACCAATAA"));
	assert(!filter2.contains("ATCGGGTCATCAACCAATTT"));
	cout << "premade bf tests done" << endl;

	//memory leak tests
	BloomFilter* filter3 = new BloomFilter(filterSize, hashMan);
	cout << memory_usage() - memUsage << "kb" << endl;
	delete(filter3);
	cout << memory_usage() - memUsage << "kb" << endl;

//	vector<vector<int> > test(1, vector<int>(1, 1));

//	vector<BloomFilter> test;
//	test.push_back(filter2);

	cout << "memory leak prevention tests done" << endl;
	cout << memory_usage() - memUsage << "kb" << endl;

	cout << "done" << endl;
	return 0;
}

