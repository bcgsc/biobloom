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

int main(int argc, char **argv)
{

	string infoFile = "Test.txt";

	vector<string> map;

	map.push_back("/home/pubseq/genomes/9606/hg19/bwa_ind/genome/human.fasta");

	BloomFilterInfo info("HG19_chr21", 33, 0.02, 47000000, map, 6);

//	//test getting Optimal Number of hash functions' function.
//	assert(info.calcOptiHashNum(16,1) == 11);

	info.printInfoFile(infoFile);

	cout << "Output tests done. check info file: " << infoFile << endl;
	while (1)
	{
	    if (getchar())
	       break;
	}

	//test loading of info into new object;
	BloomFilterInfo info2(infoFile);

	//should be fairly identical
	assert(info2.getCalcuatedFilterSize() == info.getCalcuatedFilterSize());
	assert(info2.getFilterID() == info.getFilterID());

	cout << "Asserts complete. cleaning up" << endl;
	remove(infoFile.c_str());

	return 0;
}
