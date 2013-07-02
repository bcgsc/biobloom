/*
 * WindowedParser Unit tests
 * BloomFilterGenerator Unit tests
 */

#include "BioBloomMaker/WindowedFileParser.h"
#include "BioBloomMaker/WindowedFileParser.cpp"
#include "BioBloomMaker/BloomFilterGenerator.h"
#include "BioBloomMaker/BloomFilterGenerator.cpp"
#include "Common/BloomFilter.h"
#include "Common/HashManager.h"
#include <string>
#include <assert.h>
#include <iostream>
#include <boost/unordered/unordered_map.hpp>
#include <vector>

using namespace std;

int main(int argc, char **argv) {

	//Using human genome reference as test data
	string fileName =
			"/projects/cjustin/BioBloomTools/reference/NC_000913.fasta";

	//test BloomFilterGenerator

	vector<string> filenames;
	filenames.push_back(fileName);
	//test count
	BloomFilterGenerator gen(filenames, 20);

	cout << gen.getExpectedEntries() << "bp" << endl;

	//set filter parameters
//	size_t filterSize = 34359738368;
	size_t filterSize = gen.getExpectedEntries()*8 + (64 - ((gen.getExpectedEntries()*8)% 64));

	gen.setFilterSize(filterSize);
	gen.addHashFuncs(6);
	string filename = "test.bf";

	gen.generate(filename);
	//Check storage can occur properly

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

	HashManager test;
	test.addHashFunction("CityHash64", 0);
	test.addHashFunction("CityHash64", 1);
	test.addHashFunction("CityHash64", 2);
	test.addHashFunction("CityHash64", 3);
	test.addHashFunction("CityHash64", 4);
	test.addHashFunction("CityHash64", 5);

	//check loading of stored filter
	BloomFilter filter2(filterSize, filename, test);

	//Check if loaded filter is able to report expected results
	assert(filter2.contains("AGCTTTTCATTCTGACTGCA"));
	assert(!filter2.contains("GGCTTTTCATTCTGACTGCA"));

	cout << "BloomFilterGenerator Tests Done. Cleaning up" << endl;
	remove(filename.c_str());

	return 0;
}
