/*
 * WindowedParser Unit tests
 * BloomFilterGenerator Unit tests
 */

#include "BioBloomMaker/WindowedFileParser.h"
#include "BioBloomMaker/WindowedFileParser.cpp"
#include "BioBloomMaker/BloomFilterGenerator.h"
#include "BioBloomMaker/BloomFilterGenerator.cpp"
#include "Common/BloomFilter.h"
#include <string>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <boost/unordered/unordered_map.hpp>
#include <vector>
#include <fstream>
#include "Common/ReadsProcessor.h"

using namespace std;

int main(int argc, char **argv) {

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
	} else {
		assert(filterSize == fileSize * 8);
	}
	ifile.close();

	//check loading of stored filter
	BloomFilter filter2(filterSize, hashNum, kmerSize, filename);

	ReadsProcessor proc(kmerSize);

//	//Check if loaded filter is able to report expected results

//	int countsHit = 0;
//
//	for (int i = 0; i < 100000; ++i) {
//		if(filter2.contains(proc.prepSeq("AGCTTTTCATTCTGACTGCA", 0)))
//		{
//			++countsHit;
//		}
//	}
//	cout << countsHit << endl;
	assert(filter2.contains(proc.prepSeq("AGCTTTTCATTCTGACTGCA", 0)));
	assert(!filter2.contains(proc.prepSeq("GGCTTTTCATTCTGACTGCA", 0)));

	cout << "BloomFilterGenerator Tests Done. Cleaning up" << endl;
	remove(filename.c_str());

	return 0;
}
