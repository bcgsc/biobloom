/*
 * WindowedParser Unit tests
 * BloomFilterGenerator Unit tests
 */

#include "Common/BloomFilterInfo.h"
#include <assert.h>
#include <iostream>
#include <boost/unordered/unordered_map.hpp>
#include <vector>
#include <string>
#include <cassert>

using namespace std;

int main(int argc, char **argv) {

	string infoFile = "/home/cjustin/workspace/TestData/20120820_HG19_Chr21.txt";

	vector<string> map;
	map[0] = "/home/pubseq/genomes/9606/hg19/bwa_ind/genome/human.fasta";
	BloomFilterInfo info("HG19_chr21", 33, 0.02, 47000000, map, 6);
	info.addHashFunction("CityHash64", 0);
	info.addHashFunction("CityHash64", 1);
	info.addHashFunction("CityHash64", 2);
	info.addHashFunction("CityHash64", 3);
	info.addHashFunction("CityHash64", 4);
	info.addHashFunction("CityHash64", 5);

//	//test getting Optimal Number of hash functions' function.
//	assert(info.calcOptiHashNum(16,1) == 11);

	info.printInfoFile(infoFile);

	cout << "Output tests done. check info file: " << infoFile << endl;

	//test loading of info into new object;
	BloomFilterInfo info2(infoFile);

	//should be fairly identical
	assert(info2.getCalcuatedFilterSize() == info.getCalcuatedFilterSize());
	assert(info2.getFilterID() == info.getFilterID());
	assert(info2.getSeedHashSigniture() == info.getSeedHashSigniture());

	return 0;
}
