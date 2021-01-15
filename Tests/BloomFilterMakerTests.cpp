/*
 * WindowedParser Unit tests
 * BloomFilterGenerator Unit tests
 */

#include "BioBloomMaker/BloomFilterGenerator.h"
#include "BioBloomMaker/BloomFilterGenerator.cpp"
#include <string>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <fstream>
#include "btl_bloomfilter/BloomFilter.hpp"
#include "btl_bloomfilter/vendor/ntHashIterator.hpp"
#include "omp.h"

using namespace std;

int main() {

	//Load some testdata
	string fileName = "ecoli.fasta";

	size_t kmerSize = 20;
	size_t hashNum = 5;

	//test BloomFilterGenerator

	vector<string> filenames;
	filenames.push_back(fileName);
	//test count
	BloomFilterGenerator filterGen(filenames, kmerSize, hashNum);

	cout << filterGen.getExpectedEntries() << "bp" << endl;

	//set filter parameters
//	size_t filterSize = 34359738368;
	size_t filterSize = filterGen.getExpectedEntries() * 8
			+ (64 - ((filterGen.getExpectedEntries() * 8) % 64));
	string filename = "test.bf";

	//set filter size
	filterGen.setFilterSize(filterSize);

	filterGen.generate(filename);
	//Check storage can occur properly

	ifstream ifile(filename.c_str());
	assert(ifile.is_open());
	ifile.seekg(0, ios::end); // move to end of file
	size_t fileSize = ifile.tellg(); // file size in bytes
	//file size should be same as filter size (Round to block size)
	if (filterSize % 64 > 0) {
		assert((filterSize + (64 - (filterSize % 64))) == fileSize * 8);
	}
//	else {
//		assert(filterSize == fileSize * 8);
//	}
	ifile.close();

	//check loading of stored filter
	BloomFilter filter2(filename);

	assert(filter2.contains(*ntHashIterator("AGCTTTTCATTCTGACTGCA", filter2.getHashNum(), filter2.getKmerSize())));
	assert(!filter2.contains(*ntHashIterator("AGCTTTTCATTCTGACTGCG", filter2.getHashNum(), filter2.getKmerSize())));

	cout << "BloomFilterGenerator Tests Done. Cleaning up" << endl;
	remove(filename.c_str());

	return 0;
}
