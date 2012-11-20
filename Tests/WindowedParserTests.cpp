/*
 * WindowedParser Unit tests
 * BloomFilterGenerator Unit tests
 */

#include "BioBloomMaker/WindowedFileParser.h"
#include "BioBloomMaker/WindowedFileParser.cpp"
#include "BioBloomMaker/BloomFilterGenerator.h"
#include "BioBloomMaker/BloomFilterGenerator.cpp"
#include <string>
#include <assert.h>
#include <iostream>
#include <boost/unordered/unordered_map.hpp>
#include <vector>

using namespace std;

int main(int argc, char **argv) {

	//test windowed parser
	int16_t windowSize = 50;

	//Using human genome reference as test data
	string fileName =
			"/home/pubseq/genomes/9606/hg19/bwa_ind/genome/human.fasta";
	WindowedFileParser testParser(fileName, windowSize);
	assert(testParser.getHeaders().size() == 25);

	//test setting location
	testParser.setLocationByHeader("chr1");
	testParser.setLocationByHeader("chr2");
	testParser.setLocationByHeader("chr9");
	testParser.setLocationByHeader("chrM");
	testParser.setLocationByHeader("chr3");
	testParser.setLocationByHeader("chr21");

//todo: more tests needed for windowed parser (take a look at header file)
//	string seq = testParser.getNext();
//	while (!seq.empty()) {
//		cout << seq;
//		seq = testParser.getNext();
//	}

	//Windowed Parser Tests done
	cout << "Windowed Parser Tests Done." << endl;
	return 0;
}
