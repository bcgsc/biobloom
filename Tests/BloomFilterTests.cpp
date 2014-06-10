/*
 * BloomFilterTests.cpp
 * Unit Tests for hashmanager and bloomfilter classes
 *  Created on: Aug 14, 2012
 *      Author: cjustin
 */

#include "Common/BloomFilter.h"
#include <string>
#include <assert.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include "Common/ReadsProcessor.h"
#if _OPENMP
# include <omp.h>
#endif

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
	//memory usage from before
	int memUsage = memory_usage();

	size_t filterSize = 1000000000;
	BloomFilter filter(filterSize, 5, 20);
	ReadsProcessor proc(20);

	filter.insert(proc.prepSeq("ATCGGGTCATCAACCAATAT", 0));
	filter.insert(proc.prepSeq("ATCGGGTCATCAACCAATAC", 0));
	filter.insert(proc.prepSeq("ATCGGGTCATCAACCAATAG", 0));
	filter.insert(proc.prepSeq("ATCGGGTCATCAACCAATAA", 0));

	//Check if filter is able to report expected results
	assert(filter.contains(proc.prepSeq("ATCGGGTCATCAACCAATAT", 0)));
	assert(filter.contains(proc.prepSeq("ATCGGGTCATCAACCAATAC", 0)));
	assert(filter.contains(proc.prepSeq("ATCGGGTCATCAACCAATAG", 0)));
	assert(filter.contains(proc.prepSeq("ATCGGGTCATCAACCAATAA", 0)));

	assert(!filter.contains(proc.prepSeq("ATCGGGTCATCAACCAATTA",0)));
	assert(!filter.contains(proc.prepSeq("ATCGGGTCATCAACCAATTC",0)));

	//should be size of bf (amortized)
	cout << memory_usage() - memUsage << "kb" << endl;

	cout << "de novo bf tests done" << endl;

	//Check storage can occur properly
	string filename = "/tmp/bloomFilter.bf";
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
	BloomFilter filter2(filterSize, 5, 20, filename);

	//should be double size of bf (amortized)
	cout << memory_usage() - memUsage << "kb" << endl;

	//Check if loaded filter is able to report expected results
	assert(filter2.contains(proc.prepSeq("ATCGGGTCATCAACCAATAT", 0)));
	assert(filter2.contains(proc.prepSeq("ATCGGGTCATCAACCAATAC", 0)));
	assert(filter2.contains(proc.prepSeq("ATCGGGTCATCAACCAATAG", 0)));
	assert(filter2.contains(proc.prepSeq("ATCGGGTCATCAACCAATAA", 0)));

	assert(!filter2.contains(proc.prepSeq("ATCGGGTCATCAACCAATTA",0)));
	assert(!filter2.contains(proc.prepSeq("ATCGGGTCATCAACCAATTC",0)));
	cout << "premade bf tests done" << endl;

	//memory leak tests
	BloomFilter* filter3 = new BloomFilter(filterSize, 5, 20);

	size_t tempMem = memory_usage() - memUsage;

	cout << memory_usage() - memUsage << "kb" << endl;
	delete(filter3);
	cout << memory_usage() - memUsage << "kb" << endl;
	assert(tempMem != memory_usage() - memUsage);

//	vector<vector<int> > test(1, vector<int>(1, 1));

//	vector<BloomFilter> test;
//	test.push_back(filter2);

	cout << "memory leak prevention tests done" << endl;
	cout << memory_usage() - memUsage << "kb" << endl;

	remove(filename.c_str());

//	//check parallelized code speed
//	cout << "testing code parallelization" << endl;
//	double start_s = omp_get_wtime();
//
//	cout << "parallelization" << endl;
//	for( int i =0; i < 10000000; i++)
//	{
//		vector<size_t> values = multiHash("ATCGGGTCATCAACCAATAA", 5, 20);
//		filter.contains(values);
//	}
//	double stop_s=omp_get_wtime();
//	cout << "time: " << stop_s-start_s << endl;
//
//	start_s = omp_get_wtime();
//	cout << "non parallelization" << endl;
//	for( int i =0; i < 10000000; i++)
//	{
//		vector<size_t> values = multiHashNonPara("ATCGGGTCATCAACCAATAA", 5, 20);
//		filter.contains(values);
//	}
//	stop_s=omp_get_wtime();
//	cout << "time: " << stop_s-start_s << endl;

	cout << "done" << endl;
	return 0;
}

