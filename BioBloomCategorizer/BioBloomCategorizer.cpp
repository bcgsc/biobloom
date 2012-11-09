/*
 * BioBloomCategorizer.cpp
 *
 *  Created on: Sep 7, 2012
 *      Author: cjustin
 */
//todo: refactor into additional class
//todo: UNIT TESTS!
#include <sstream>
#include <string>
#include <getopt.h>
#include "boost/unordered/unordered_map.hpp"
#include "Common/BloomFilterInfo.h"
#include "Common/HashManager.h"
#include "MultiFilter.h"
#include <vector>
#include "DataLayer/FastaReader.h"
#include <iostream>
#include <fstream>
#include "boost/shared_ptr.hpp"
#include "Common/ReadsProcessor.h"
#include "Common/Uncompress.h"
#include <sys/stat.h>
#include <sys/types.h>
using namespace std;

//group filters with same hash signature
//todo: use with auto_ptr or unique_ptr?
boost::unordered_map<string, vector<boost::shared_ptr<BloomFilterInfo> > > infoFiles;
boost::unordered_map<string, boost::shared_ptr<MultiFilter> > filters;
vector<string> hashSigs;

/*
 * Parses input string into seperate strings, returning a vector.
 */
vector<string> convertInputString(const string &inputString) {
	vector<string> currentInfoFile;
	string temp;
	stringstream converter(inputString);
	while (converter >> temp) {
		currentInfoFile.push_back(temp);
	}
	return currentInfoFile;
}

/*
 * checks if file exists
 */
bool fexists(const string &filename) {
	ifstream ifile(filename.c_str());
	return ifile;
}

void loadFilters(const vector<string> &filterFilePaths) {
	cout << "Starting to Load Filters." << endl;
	//load up files
	for (vector<string>::const_iterator it = filterFilePaths.begin();
			it != filterFilePaths.end(); ++it)
	{
		//check if files exist
		if (!fexists(*it)) {
			cerr << "Error: " + (*it) + " File cannot be opened" << endl;
			exit(1);
		}
		string infoFileName = (*it).substr(0, (*it).length() - 2) + "txt";
		if (!fexists(infoFileName)) {
			cerr
					<< "Error: " + (infoFileName)
							+ " File cannot be opened. A corresponding info file is needed."
					<< endl;
			exit(1);
		}

		//info file creation
		boost::shared_ptr<BloomFilterInfo> info(
				new BloomFilterInfo(infoFileName));
		//append kmer size to hash signature to insure correct kmer size is used
		stringstream hashSig;
		hashSig << info->getKmerSize() << info->getSeedHashSigniture();

		//if hashSig exists add filter to list
		if (infoFiles.count(hashSig.str()) == 1) {
			infoFiles[hashSig.str()].push_back(info);
			filters[hashSig.str()]->addFilter(info->getFilterID(), *it);
		} else {
			vector<boost::shared_ptr<BloomFilterInfo> > tempVect;
			tempVect.push_back(info);
			vector<size_t>::const_iterator seeds =
					info->getSeedValues().begin();
			//Create HashManager for MultiFilter
			HashManager hashMan;
			for (vector<string>::const_iterator hashFn =
					info->getHashFunctionNames().begin();
					hashFn != info->getHashFunctionNames().end(); ++hashFn)
			{
				hashMan.addHashFunction(*hashFn, *seeds);
				++seeds;
			}
			hashSigs.push_back(hashSig.str());
			boost::shared_ptr<MultiFilter> temp(new MultiFilter(hashMan));
			filters[hashSig.str()] = temp;
			filters[hashSig.str()]->addFilter(info->getFilterID(), *it);
			infoFiles[hashSig.str()] = tempVect;
		}
		cout << "Loaded Filter: " + info->getFilterID() << endl;
	}
	cout << "Filter Loading Complete." << endl;
}

int main(int argc, char *argv[]) {

	//switch statement variable
	int c;

	//command line variables
	string rawInputFiles = "";
	string outputDir = "";
	string filtersFile = "";
	int16_t minHit = 2;
	bool printReads = false;

	//long form arguments
	//each option format { "optionName", necessary option or not, I have no idea, 'symbol'}
	static struct option long_options[] = { { "input_files", 1, NULL, 'i' }, {
			"output_dir", 0, NULL, 'o' }, { "filters", 1, NULL, 'f' }, {
			"min_hit", 0, NULL, 'm' }, { "print_fasta", 1, NULL, 'p' }, { NULL,
			0, NULL, 0 } };

	//actual checking step
	int option_index = 0;
	while ((c = getopt_long(argc, argv, "i:f:o:k:m:p", long_options,
			&option_index)) != -1)
	{
		switch (c) {
		case 'i': {
			rawInputFiles = optarg;
			break;
		}
		case 'f': {
			filtersFile = optarg;
			break;
		}
		case 'o': {
			outputDir = optarg;
			if (outputDir.at(outputDir.length() - 1) != '/') {
				outputDir = outputDir + '/';
			}
			break;
		}
		case 'm': {
			stringstream convert(optarg);
			if (!(convert >> minHit)) {
				cout << "Error - Invalid parameter! m: " << optarg << endl;
				return 0;
			}
			break;
		}
		case 'p': {
			printReads = true;
			break;
		}
		default: {
			exit(1);
		}
		}
	}

	vector<string> filterFilePaths = convertInputString(filtersFile);
	vector<string> inputFiles = convertInputString(rawInputFiles);

	loadFilters(filterFilePaths);

	//filtering step
	//create directory structure if it does not exist
	if (!fexists(outputDir)) {
		mkdir(outputDir.c_str(), 0755);
	}

	ofstream summaryOutput((outputDir + "summary.tsv").c_str(), ios::out);
	ofstream readStatusOutput((outputDir + "readStatus.tsv").c_str(), ios::out);

	//print header
	summaryOutput << "type";
	readStatusOutput << "readID\tseqSize";

	//variables for storing results summary
	boost::unordered_map<string, size_t> aboveThreshosld;
	boost::unordered_map<string, size_t> belowThreshosld;
	size_t totalReads = 0;

	//output file bins
	if (printReads) {
		boost::unordered_map<string, boost::shared_ptr<ofstream> > outputFiles;
		boost::shared_ptr<ofstream> noMatch(
				new ofstream((outputDir + "noMatch.fasta").c_str(), ios::out));
		boost::shared_ptr<ofstream> multiMatch(
				new ofstream((outputDir + "multiMatch.fasta").c_str(),
						ios::out));
		outputFiles["noMatch"] = noMatch;
		outputFiles["multiMatch"] = multiMatch;

		//initialize variables and print filter ids
		for (vector<string>::iterator j = hashSigs.begin(); j != hashSigs.end();
				++j)
		{
			vector<string> idsInFilter = (*filters[*j]).getFilterIds();
			for (vector<string>::iterator i = idsInFilter.begin();
					i != idsInFilter.end(); ++i)
			{
				boost::shared_ptr<ofstream> temp(
						new ofstream((outputDir + *i + ".fasta").c_str(),
								ios::out));
				outputFiles[*i] = temp;
				summaryOutput << "\t" << *i << "_"
						<< (*(infoFiles[*j].front())).getKmerSize();
				readStatusOutput << "\t" << *i << "_"
						<< (*(infoFiles[*j].front())).getKmerSize();
				aboveThreshosld[*i] = 0;
				belowThreshosld[*i] = 0;
			}
		}
		summaryOutput << "\n";
		readStatusOutput << "\n";

		//Todo: make sure this prints out only when filters are loaded
		//gcc currently optimizes to print this before loading can complete
		cout << "Filtering Start" << endl;

		for (vector<string>::iterator it = inputFiles.begin();
				it != inputFiles.end(); ++it)
		{
			FastaReader sequence((*it).c_str(), FastaReader::NO_FOLD_CASE);
			FastqRecord rec;
			//hits results stored in hashmap of filternames and hits
			boost::unordered_map<string, size_t> hits(filterFilePaths.size());
			while (sequence >> rec) {
				//split reads into kmerSizes specified (ignore trailing bases)

				//for skipping bad reads
				bool readOK = true;

				//initialize hits to zero
				for (vector<string>::iterator j = hashSigs.begin();
						j != hashSigs.end(); ++j)
				{
					const vector<string> &idsInFilter =
							(*filters[*j]).getFilterIds();
					for (vector<string>::const_iterator i = idsInFilter.begin();
							i != idsInFilter.end(); ++i)
					{
						hits[*i] = 0;
					}
				}

				//for each hashSigniture/kmer combo multi, cut up read into kmer sized used
				for (vector<string>::iterator j = hashSigs.begin();
						j != hashSigs.end(); ++j)
				{
					//get filterIDs to iterate through has in a consistent order
					const vector<string> &idsInFilter =
							(*filters[*j]).getFilterIds();

					//get kmersize for set of info files
					int16_t kmerSize = (*(infoFiles[*j].front())).getKmerSize();
					size_t currentKmerNum = 0;

					//Establish tiling pattern
					int16_t startModifier = (rec.seq.length() % kmerSize) / 2;

					ReadsProcessor proc(kmerSize);
					//cut read into kmer size given
					while (rec.seq.length() >= (currentKmerNum + 1) * kmerSize)
					{

						const string &currentKmer = proc.prepSeq(rec.seq,
								currentKmerNum * kmerSize + startModifier);

						//check to see if string is invalid
						if (!currentKmer.empty()) {
							const boost::unordered_map<string, bool> &results =
									filters[*j]->multiContains(currentKmer);

							//record hit number in order
							for (vector<string>::const_iterator i =
									idsInFilter.begin(); i != idsInFilter.end();
									++i)
							{
								if (results.find(*i)->second) {
									++hits[*i];
								}
							}
						} else {
							readOK = false;
							break;
						}
						++currentKmerNum;
					}
				}

				if (readOK) {
					//print readID
					readStatusOutput << rec.id << "\t" << rec.seq.length();
					++totalReads;
					if (totalReads % 1000000 == 0) {
						cout << "Currently Reading Read Number: " << totalReads
								<< endl;
					}
					int16_t totalHits = 0;
					for (vector<string>::iterator j = hashSigs.begin();
							j != hashSigs.end(); ++j)
					{
						//update summary
						const vector<string> &idsInFilter =
								(*filters[*j]).getFilterIds();
						for (vector<string>::const_iterator i =
								idsInFilter.begin(); i != idsInFilter.end();
								++i)
						{
							readStatusOutput << "\t" << hits[*i];
							if (hits[*i] >= minHit) {
								++totalHits;
								++aboveThreshosld[*i];
							} else if (hits[*i] != 0) {
								++belowThreshosld[*i];
							}
						}
					}
					if (totalHits == 0) {
						(*outputFiles["noMatch"]) << ">" << rec.id << "\n"
								<< rec.seq << endl;
					} else if (totalHits > 1) {
						(*outputFiles["multiMatch"]) << ">" << rec.id << "\n"
								<< rec.seq << endl;
					} else {
						for (vector<string>::iterator j = hashSigs.begin();
								j != hashSigs.end(); ++j)
						{
							vector<string> idsInFilter =
									(*filters[*j]).getFilterIds();
							for (vector<string>::iterator i =
									idsInFilter.begin(); i != idsInFilter.end();
									++i)
							{
								if (hits[*i] >= minHit) {
									(*outputFiles[*i]) << ">" << rec.id << "\n"
											<< rec.seq << endl;
									break;
								}
							}
						}
					}
					readStatusOutput << "\n";
				}
			}
		}

		//close sorting files
		for (vector<string>::iterator j = hashSigs.begin(); j != hashSigs.end();
				++j)
		{
			vector<string> idsInFilter = (*filters[*j]).getFilterIds();
			for (vector<string>::iterator i = idsInFilter.begin();
					i != idsInFilter.end(); ++i)
			{
				outputFiles[*i]->close();
			}
		}
		outputFiles["noMatch"]->close();
		outputFiles["multiMatch"]->close();
	} else {
		//initialize variables and print filter ids
		for (vector<string>::iterator j = hashSigs.begin(); j != hashSigs.end();
				++j)
		{
			vector<string> idsInFilter = (*filters[*j]).getFilterIds();
			for (vector<string>::iterator i = idsInFilter.begin();
					i != idsInFilter.end(); ++i)
			{
				summaryOutput << "\t" << *i << "_"
						<< (*(infoFiles[*j].front())).getKmerSize();
				readStatusOutput << "\t" << *i << "_"
						<< (*(infoFiles[*j].front())).getKmerSize();
				aboveThreshosld[*i] = 0;
				belowThreshosld[*i] = 0;
			}
		}
		summaryOutput << "\n";
		readStatusOutput << "\n";

		//Todo: make sure this prints out only when filters are loaded
		//gcc currently optimizes to print this before loading can complete
		cout << "Filtering Start" << endl;

		for (vector<string>::iterator it = inputFiles.begin();
				it != inputFiles.end(); ++it)
		{
			FastaReader sequence((*it).c_str(), FastaReader::NO_FOLD_CASE);
			FastqRecord rec;
			//hits results stored in hashmap of filternames and hits
			boost::unordered_map<string, size_t> hits(filterFilePaths.size());
			while (sequence >> rec) {
				//split reads into kmerSizes specified (ignore trailing bases)

				//for skipping bad reads
				bool readOK = true;

				//initialize hits to zero
				for (vector<string>::iterator j = hashSigs.begin();
						j != hashSigs.end(); ++j)
				{
					const vector<string> &idsInFilter =
							(*filters[*j]).getFilterIds();
					for (vector<string>::const_iterator i = idsInFilter.begin();
							i != idsInFilter.end(); ++i)
					{
						hits[*i] = 0;
					}
				}

				//for each hashSigniture/kmer combo multi, cut up read into kmer sized used
				for (vector<string>::iterator j = hashSigs.begin();
						j != hashSigs.end(); ++j)
				{
					//get filterIDs to iterate through has in a consistent order
					const vector<string> &idsInFilter =
							(*filters[*j]).getFilterIds();

					//get kmersize for set of info files
					int16_t kmerSize = (*(infoFiles[*j].front())).getKmerSize();
					size_t currentKmerNum = 0;

					//Establish tiling pattern
					int16_t startModifier = (rec.seq.length() % kmerSize) / 2;

					ReadsProcessor proc(kmerSize);
					//cut read into kmer size given
					while (rec.seq.length() >= (currentKmerNum + 1) * kmerSize)
					{

						const string &currentKmer = proc.prepSeq(rec.seq,
								currentKmerNum * kmerSize + startModifier);

						//check to see if string is invalid
						if (!currentKmer.empty()) {
							const boost::unordered_map<string, bool> &results =
									filters[*j]->multiContains(currentKmer);

							//record hit number in order
							for (vector<string>::const_iterator i =
									idsInFilter.begin(); i != idsInFilter.end();
									++i)
							{
								if (results.find(*i)->second) {
									++hits[*i];
								}
							}
						} else {
							readOK = false;
							break;
						}
						++currentKmerNum;
					}
				}

				if (readOK) {
					//print readID
					readStatusOutput << rec.id << "\t" << rec.seq.length();
					++totalReads;
					if (totalReads % 1000000 == 0) {
						cout << "Currently Reading Read Number: " << totalReads
								<< endl;
					}
					int16_t totalHits = 0;
					for (vector<string>::iterator j = hashSigs.begin();
							j != hashSigs.end(); ++j)
					{
						//update summary
						const vector<string> &idsInFilter =
								(*filters[*j]).getFilterIds();
						for (vector<string>::const_iterator i =
								idsInFilter.begin(); i != idsInFilter.end();
								++i)
						{
							readStatusOutput << "\t" << hits[*i];
							if (hits[*i] >= minHit) {
								++totalHits;
								++aboveThreshosld[*i];
							} else if (hits[*i] != 0) {
								++belowThreshosld[*i];
							}
						}
					}
					readStatusOutput << "\n";
				}
			}
		}
	}
	readStatusOutput << endl;
	readStatusOutput.close();

	//print summary information and close filehandles
	summaryOutput << ">=" << minHit << "_proportion";
	for (vector<string>::iterator j = hashSigs.begin(); j != hashSigs.end();
			++j)
	{
		vector<string> idsInFilter = (*filters[*j]).getFilterIds();
		for (vector<string>::iterator i = idsInFilter.begin();
				i != idsInFilter.end(); ++i)
		{
			summaryOutput << "\t"
					<< double(aboveThreshosld[*i]) / double(totalReads);
		}
	}
	summaryOutput << "\n<" << minHit << "_proportion";
	for (vector<string>::iterator j = hashSigs.begin(); j != hashSigs.end();
			++j)
	{
		vector<string> idsInFilter = (*filters[*j]).getFilterIds();
		for (vector<string>::iterator i = idsInFilter.begin();
				i != idsInFilter.end(); ++i)
		{
			summaryOutput << "\t"
					<< double(belowThreshosld[*i]) / double(totalReads);
		}
	}
	summaryOutput << "\n" << "0_proportion";
	for (vector<string>::iterator j = hashSigs.begin(); j != hashSigs.end();
			++j)
	{
		vector<string> idsInFilter = (*filters[*j]).getFilterIds();
		for (vector<string>::iterator i = idsInFilter.begin();
				i != idsInFilter.end(); ++i)
		{
			summaryOutput << "\t"
					<< double(
							totalReads - belowThreshosld[*i]
									- aboveThreshosld[*i]) / double(totalReads);
		}
	}

	summaryOutput << "\n>=" << minHit << "_reads";
	for (vector<string>::iterator j = hashSigs.begin(); j != hashSigs.end();
			++j)
	{
		vector<string> idsInFilter = (*filters[*j]).getFilterIds();
		for (vector<string>::iterator i = idsInFilter.begin();
				i != idsInFilter.end(); ++i)
		{
			summaryOutput << "\t" << aboveThreshosld[*i];
		}
	}
	summaryOutput << "\n<" << minHit << "_reads";
	for (vector<string>::iterator j = hashSigs.begin(); j != hashSigs.end();
			++j)
	{
		vector<string> idsInFilter = (*filters[*j]).getFilterIds();
		for (vector<string>::iterator i = idsInFilter.begin();
				i != idsInFilter.end(); ++i)
		{
			summaryOutput << "\t" << belowThreshosld[*i];
		}
	}
	summaryOutput << "\n" << "0_reads";
	for (vector<string>::iterator j = hashSigs.begin(); j != hashSigs.end();
			++j)
	{
		vector<string> idsInFilter = (*filters[*j]).getFilterIds();
		for (vector<string>::iterator i = idsInFilter.begin();
				i != idsInFilter.end(); ++i)
		{
			summaryOutput << "\t"
					<< totalReads - belowThreshosld[*i] - aboveThreshosld[*i];
		}
	}
	summaryOutput.close();
	cout << "Total Reads:" << totalReads << endl;
}
