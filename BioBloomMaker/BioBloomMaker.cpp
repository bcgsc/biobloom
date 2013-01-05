/*
 * BioBloom.cpp
 *
 *  Created on: Jun 20, 2012
 *      Author: cjustin
 */

#include <sstream>
#include <string>
#include <vector>
#include "BloomFilterGenerator.h"
#include "Common/BloomFilterInfo.h"
#include <boost/unordered/unordered_map.hpp>
#include <getopt.h>
using namespace std;

void printHelpDialog()
{
	static const char dialog[] =
			"Usage: BioBloomMaker -p [FILTERID] [OPTION]... [FILE]...\n"
					"Creates a bf and txt file from a list of fasta files. The input sequences are\n"
					"cut into a k-mers with a sliding window and their hash signatures are inserted\n"
					"into a bloom filter.\n"
					"\n"
					"  -f, --fal_pos_rate=N   Maximum false positive rate to use in filter. [0.02]\n"
					"  -p, --file_prefix=N    Filter prefix and filter ID. Required option.\n"
					"  -o, --output_dir=N     Output location of the filter and filter info files.\n"
					"  -g, --hash_num=N       Set number of hash functions to use in filter instead\n"
					"                         of automatically using calculated optimal number of\n"
					"                         functions.\n"
					"  -k, --output_fastq     K-mer size to use to create filter.[25]\n"
					"  -h, --help             Display this dialog.\n"
					"\n"
					"Report bugs to <cjustin@bcgsc.ca>.";
	cerr << dialog << endl;
	exit(0);
}

int main(int argc, char *argv[])
{

	bool die = false;

	//switch statement variable
	int c;

	//command line variables
	double fpr = 0.02;
	string filterPrefix = "";
	string outputDir = "";
	int16_t kmerSize = 25;
	int16_t hashNum = 0;

	//long form arguments
	//each option format { "optionName", necessary option or not, I have no idea, 'symbol'}
	static struct option long_options[] = {
			{
					"fal_pos_rate", required_argument, NULL, 'f' }, {
					"file_prefix", required_argument, NULL, 'p' }, {
					"output_dir", optional_argument, NULL, 'o' }, {
					"hash_num", 0, NULL, 'g' }, {
					"kmer_size", 1, NULL, 'k' }, {
					"help", no_argument, NULL, 'h' }, {
					NULL, 0, NULL, 0 } };

	//actual checking step
	int option_index = 0;
	while ((c = getopt_long(argc, argv, "f:p:o:k:n:g:h", long_options,
			&option_index)) != -1)
	{
		switch (c) {
		case 'f': {
			stringstream convert(optarg);
			if (!(convert >> fpr)) {
				cerr << "Error - Invalid set of bloom filter parameters! f: "
						<< optarg << endl;
				return 0;
			}
			if (fpr > 1) {
				cerr << "Error -f cannot be greater than 1 " << optarg << endl;
				return 0;
			}
			break;
		}
		case 'p': {
			filterPrefix = optarg;
			break;
		}
		case 'o': {
			outputDir = optarg;
			if (outputDir.at(outputDir.length() - 1) != '/') {
				outputDir = outputDir + '/';
			}
			break;
		}
		case 'k': {
			stringstream convert(optarg);
			if (!(convert >> kmerSize)) {
				cerr << "Error - Invalid set of bloom filter parameters! k: "
						<< optarg << endl;
				return 0;
			}
			break;
		}
		case 'g': {
			stringstream convert(optarg);
			if (!(convert >> hashNum)) {
				cerr << "Error - Invalid set of bloom filter parameters! g: "
						<< optarg << endl;
				return 0;
			}
			break;
		}
		case 'h': {
			printHelpDialog();
			break;
		}
		default: {
			die = true;
			break;
		}
		}
	}

	//Stores fasta input file names
	vector<string> inputFiles;

	while (optind < argc) {
		inputFiles.push_back(argv[optind]);
		optind++;
	}

	//Check needed options
	if (inputFiles.size() == 0) {
		cerr << "Need Input File" << endl;
		die = true;
	}
	if (filterPrefix.size() == 0) {
		cerr << "Need Filter Prefix ID" << endl;
		die = true;
	}
	if (filterPrefix.find('/') != string::npos) {
		cerr << "Prefix ID cannot have '/' characters" << endl;
		die = true;
	}
	if (die) {
		cerr << "Try '--help' for more information.\n";
		exit(EXIT_FAILURE);
	}

	//create filter
	BloomFilterGenerator filterGen(inputFiles, kmerSize);

	//get expected number of entries
	size_t entryNum = filterGen.getExpectedEntries();

	//set number of hash functions used
	if (hashNum == 0) {
		//get optimal number of hash functions
		hashNum = filterGen.calcOptiHashNum(fpr);
	}

	BloomFilterInfo info(filterPrefix, kmerSize, fpr, entryNum, inputFiles,
			hashNum);

	//get calculated size of Filter
	size_t filterSize = info.getCalcuatedFilterSize();

	//Add seed to bloom filter and get them for info file
	vector<size_t> seeds = filterGen.addHashFuncs(hashNum);
	assert(seeds.size() == hashNum);
	filterGen.setFilterSize(filterSize);

	//output filter
	filterGen.generate(outputDir + filterPrefix + ".bf");
	//set hash function info
	vector<string> hashFunctions = filterGen.getHashFuncNames();

	for (int16_t i = 0; i < hashNum; ++i) {
		info.addHashFunction(hashFunctions[i], seeds[i]);
	}

	//output info
	info.printInfoFile(outputDir + filterPrefix + ".txt");
	cerr << "Filter Creation Complete" << endl;

	return 0;
}
