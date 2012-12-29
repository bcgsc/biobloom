/*
 * BioBloomClassifier.cpp
 *
 *  Created on: Oct 17, 2012
 *      Author: cjustin
 */

#include "BioBloomClassifier.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/stat.h>

BioBloomClassifier::BioBloomClassifier(const vector<string> &filterFilePaths,
		int16_t minHit, double percentMinHit) :
		minHit(minHit), percentMinHit(percentMinHit), filterNum(
				filterFilePaths.size())
{
	loadFilters(filterFilePaths);
}

void BioBloomClassifier::loadFilters(const vector<string> &filterFilePaths)
{
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
			filters[hashSig.str()]->addFilter(info->getCalcuatedFilterSize(),
					info->getFilterID(), *it);
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
			filters[hashSig.str()]->addFilter(info->getCalcuatedFilterSize(),
					info->getFilterID(), *it);
			infoFiles[hashSig.str()] = tempVect;
		}
		cerr << "Loaded Filter: " + info->getFilterID() << endl;
	}
	cerr << "Filter Loading Complete." << endl;
}

void BioBloomClassifier::filter(const vector<string> &inputFiles,
		const string &outputPrefix)
{

	ofstream readStatusOutput((outputPrefix + "_status.tsv").c_str(), ios::out);

	//print header
	readStatusOutput << "readID\tseqSize";

	//variables for storing results summary
	boost::unordered_map<string, size_t> aboveThreshold;
	boost::unordered_map<string, size_t> belowThreshold;
	size_t totalReads = 0;

	//initialize variables and print filter ids
	for (vector<string>::const_iterator j = hashSigs.begin();
			j != hashSigs.end(); ++j)
	{
		const vector<string> idsInFilter = (*filters[*j]).getFilterIds();
		for (vector<string>::const_iterator i = idsInFilter.begin();
				i != idsInFilter.end(); ++i)
		{
			readStatusOutput << "\t" << *i << "_"
					<< (*(infoFiles[*j].front())).getKmerSize();
			aboveThreshold[*i] = 0;
			belowThreshold[*i] = 0;
		}
	}
	readStatusOutput << "\n";

	//Todo: make sure this prints out only when filters are loaded
	//gcc currently optimizes to print this before loading can complete
	cerr << "Filtering Start" << endl;

	for (vector<string>::const_iterator it = inputFiles.begin();
			it != inputFiles.end(); ++it)
	{
		FastaReader sequence((*it).c_str(), FastaReader::NO_FOLD_CASE);
		FastqRecord rec;
		//hits results stored in hashmap of filternames and hits
		boost::unordered_map<string, size_t> hits(filterNum);
		while (sequence >> rec) {
			//split reads into kmerSizes specified (ignore trailing bases)

			//for skipping bad reads
			bool readOK = true;

			//initialize hits to zero
			for (vector<string>::const_iterator j = hashSigs.begin();
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
			for (vector<string>::const_iterator j = hashSigs.begin();
					j != hashSigs.end(); ++j)
			{
				if (!evaluateRead(rec, *j, hits)) {
					readOK = false;
					break;
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
				for (vector<string>::const_iterator j = hashSigs.begin();
						j != hashSigs.end(); ++j)
				{
					//update summary
					const vector<string> &idsInFilter =
							(*filters[*j]).getFilterIds();
					for (vector<string>::const_iterator i = idsInFilter.begin();
							i != idsInFilter.end(); ++i)
					{
						//print to file
						readStatusOutput << "\t" << hits[*i];

						//pick threshold, by percent or by absolute value
						int16_t kmerSize =
								(*(infoFiles[*j].front())).getKmerSize();
						size_t threshold = size_t(
								percentMinHit * (rec.seq.length() / kmerSize));
						if (minHit > threshold) {
							threshold = minHit;
						}

						if (hits[*i] >= threshold) {
							++totalHits;
							++aboveThreshold[*i];
						} else if (hits[*i] != 0) {
							++belowThreshold[*i];
						}
					}
				}
				readStatusOutput << "\n";
			}
		}
	}
	cout << "Total Reads:" << totalReads << endl;
	printSummary(outputPrefix, aboveThreshold, belowThreshold, totalReads);
}

void BioBloomClassifier::filterPrintReads(const vector<string> &inputFiles,
		const string &outputPrefix)
{

	ofstream readStatusOutput((outputPrefix + "_status.tsv").c_str(), ios::out);

	//print header
	readStatusOutput << "readID\tseqSize";

	//variables for storing results summary
	boost::unordered_map<string, size_t> aboveThreshold;
	boost::unordered_map<string, size_t> belowThreshold;
	size_t totalReads = 0;

	boost::unordered_map<string, boost::shared_ptr<ofstream> > outputFiles;
	boost::shared_ptr<ofstream> noMatch(
			new ofstream((outputPrefix + "noMatch.fastq").c_str(), ios::out));
	boost::shared_ptr<ofstream> multiMatch(
			new ofstream((outputPrefix + "multiMatch.fastq").c_str(),
					ios::out));
	outputFiles["noMatch"] = noMatch;
	outputFiles["multiMatch"] = multiMatch;

	//initialize variables and print filter ids
	for (vector<string>::const_iterator j = hashSigs.begin();
			j != hashSigs.end(); ++j)
	{
		const vector<string> idsInFilter = (*filters[*j]).getFilterIds();
		for (vector<string>::const_iterator i = idsInFilter.begin();
				i != idsInFilter.end(); ++i)
		{
			boost::shared_ptr<ofstream> temp(
					new ofstream((outputPrefix + *i + ".fastq").c_str(),
							ios::out));
			outputFiles[*i] = temp;
			readStatusOutput << "\t" << *i << "_"
					<< (*(infoFiles[*j].front())).getKmerSize();
			aboveThreshold[*i] = 0;
			belowThreshold[*i] = 0;
		}
	}
	readStatusOutput << "\n";

	//Todo: make sure this prints out only when filters are loaded
	//gcc currently optimizes to print this before loading can complete
	cerr << "Filtering Start" << endl;

	for (vector<string>::const_iterator it = inputFiles.begin();
			it != inputFiles.end(); ++it)
	{
		FastaReader sequence((*it).c_str(), FastaReader::NO_FOLD_CASE);
		FastqRecord rec;
		//hits results stored in hashmap of filternames and hits
		boost::unordered_map<string, size_t> hits(filterNum);
		while (sequence >> rec) {
			//split reads into kmerSizes specified (ignore trailing bases)

			//for skipping bad reads
			bool readOK = true;

			//initialize hits to zero
			for (vector<string>::const_iterator j = hashSigs.begin();
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
			for (vector<string>::const_iterator j = hashSigs.begin();
					j != hashSigs.end(); ++j)
			{
				if (!evaluateRead(rec, *j, hits)) {
					readOK = false;
					break;
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
				for (vector<string>::const_iterator j = hashSigs.begin();
						j != hashSigs.end(); ++j)
				{
					//update summary
					const vector<string> &idsInFilter =
							(*filters[*j]).getFilterIds();
					for (vector<string>::const_iterator i = idsInFilter.begin();
							i != idsInFilter.end(); ++i)
					{
						//print read status
						readStatusOutput << "\t" << hits[*i];

						//pick threshold, by percent or by absolute value
						int16_t kmerSize =
								(*(infoFiles[*j].front())).getKmerSize();
						size_t threshold = size_t(
								percentMinHit * (rec.seq.length() / kmerSize));
						if (minHit > threshold) {
							threshold = minHit;
						}

						if (hits[*i] >= threshold) {
							++totalHits;
							++aboveThreshold[*i];
						} else if (hits[*i] != 0) {
							++belowThreshold[*i];
						}
					}
				}
				if (totalHits == 0) {
					(*outputFiles["noMatch"]) << "@" << rec.id << "\n"
							<< rec.seq << "\n+\n" << rec.qual << endl;
				} else if (totalHits > 1) {
					(*outputFiles["multiMatch"]) << "@" << rec.id << "\n"
							<< rec.seq << "\n+\n" << rec.qual << endl;
				} else {
					for (vector<string>::const_iterator j = hashSigs.begin();
							j != hashSigs.end(); ++j)
					{
						const vector<string> idsInFilter =
								(*filters[*j]).getFilterIds();
						for (vector<string>::const_iterator i =
								idsInFilter.begin(); i != idsInFilter.end();
								++i)
						{
							//pick threshold, by percent or by absolute value
							int16_t kmerSize =
									(*(infoFiles[*j].front())).getKmerSize();
							size_t threshold = size_t(
									percentMinHit
											* (rec.seq.length() / kmerSize));
							if (minHit > threshold) {
								threshold = minHit;
							}

							if (hits[*i] >= threshold) {
								(*outputFiles[*i]) << "@" << rec.id << "\n"
										<< rec.seq << "\n+\n" << rec.qual
										<< endl;
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
	for (vector<string>::const_iterator j = hashSigs.begin();
			j != hashSigs.end(); ++j)
	{
		const vector<string> idsInFilter = (*filters[*j]).getFilterIds();
		for (vector<string>::const_iterator i = idsInFilter.begin();
				i != idsInFilter.end(); ++i)
		{
			outputFiles[*i]->close();
		}
	}
	outputFiles["noMatch"]->close();
	outputFiles["multiMatch"]->close();
	cout << "Total Reads:" << totalReads << endl;
	printSummary(outputPrefix, aboveThreshold, belowThreshold, totalReads);
}

/*
 * Prints summary information:
 * -total reads over/under threshold
 * -percent reads over threshold
 * -total reads that don't hit filter at all
 */
void BioBloomClassifier::printSummary(const string &outputPrefix,
		boost::unordered_map<string, size_t> &aboveThreshold,
		boost::unordered_map<string, size_t> &belowThreshold, size_t totalReads)
{
	ofstream summaryOutput((outputPrefix + "_summary.tsv").c_str(), ios::out);
	summaryOutput << "type";
	//initialize variables and print filter ids
	for (vector<string>::const_iterator j = hashSigs.begin();
			j != hashSigs.end(); ++j)
	{
		const vector<string> idsInFilter = (*filters[*j]).getFilterIds();
		for (vector<string>::const_iterator i = idsInFilter.begin();
				i != idsInFilter.end(); ++i)
		{
			summaryOutput << "\t" << *i << "_"
					<< (*(infoFiles[*j].front())).getKmerSize();
		}
	}
	summaryOutput << "\n";

	//print summary information and close filehandles
	summaryOutput << "\nHits";
	for (vector<string>::const_iterator j = hashSigs.begin();
			j != hashSigs.end(); ++j)
	{
		const vector<string> idsInFilter = (*filters[*j]).getFilterIds();
		for (vector<string>::const_iterator i = idsInFilter.begin();
				i != idsInFilter.end(); ++i)
		{
			summaryOutput << "\t"
					<< double(aboveThreshold[*i]) / double(totalReads);
		}
	}
	summaryOutput << "\nMiss";
	for (vector<string>::const_iterator j = hashSigs.begin();
			j != hashSigs.end(); ++j)
	{
		const vector<string> idsInFilter = (*filters[*j]).getFilterIds();
		for (vector<string>::const_iterator i = idsInFilter.begin();
				i != idsInFilter.end(); ++i)
		{
			summaryOutput << "\t"
					<< double(belowThreshold[*i]) / double(totalReads);
		}
	}
	summaryOutput << "\nConfidentMiss";
	for (vector<string>::const_iterator j = hashSigs.begin();
			j != hashSigs.end(); ++j)
	{
		const vector<string> idsInFilter = (*filters[*j]).getFilterIds();
		for (vector<string>::const_iterator i = idsInFilter.begin();
				i != idsInFilter.end(); ++i)
		{
			summaryOutput << "\t"
					<< double(
							totalReads - belowThreshold[*i]
									- aboveThreshold[*i]) / double(totalReads);
		}
	}

	summaryOutput << "\n\nHits";
	for (vector<string>::const_iterator j = hashSigs.begin();
			j != hashSigs.end(); ++j)
	{
		const vector<string> idsInFilter = (*filters[*j]).getFilterIds();
		for (vector<string>::const_iterator i = idsInFilter.begin();
				i != idsInFilter.end(); ++i)
		{
			summaryOutput << "\t" << aboveThreshold[*i];
		}
	}
	summaryOutput << "\nMiss";
	for (vector<string>::const_iterator j = hashSigs.begin();
			j != hashSigs.end(); ++j)
	{
		const vector<string> idsInFilter = (*filters[*j]).getFilterIds();
		for (vector<string>::const_iterator i = idsInFilter.begin();
				i != idsInFilter.end(); ++i)
		{
			summaryOutput << "\t" << belowThreshold[*i];
		}
	}
	summaryOutput << "\nConfidentMiss";
	for (vector<string>::const_iterator j = hashSigs.begin();
			j != hashSigs.end(); ++j)
	{
		const vector<string> idsInFilter = (*filters[*j]).getFilterIds();
		for (vector<string>::const_iterator i = idsInFilter.begin();
				i != idsInFilter.end(); ++i)
		{
			summaryOutput << "\t"
					<< totalReads - belowThreshold[*i] - aboveThreshold[*i];
		}
	}
	summaryOutput.close();

}

/*
 * Filters reads -> uses paired end information
 * Assumes only one hash signature exists (load only filters with same
 * hash functions)
 */
void BioBloomClassifier::filterPairPrint(const string &file1, const string &file2,
		const string &outputPrefix)
{
	//check if output folder exists
	string tempStr = outputPrefix.substr(0, outputPrefix.find_last_of("/"));
	folderCheck(tempStr);

	ofstream readStatusOutput((outputPrefix + "_status.tsv").c_str(), ios::out);

	//print header
	readStatusOutput << "readID\tseqSize";

	//variables for storing results summary
	boost::unordered_map<string, size_t> aboveThreshold;
	boost::unordered_map<string, size_t> belowThreshold;
	size_t totalReads = 0;

	boost::unordered_map<string, boost::shared_ptr<ofstream> > outputFiles;
	boost::shared_ptr<ofstream> noMatch1(
			new ofstream((outputPrefix + "noMatch1.fastq").c_str(), ios::out));
	boost::shared_ptr<ofstream> noMatch2(
			new ofstream((outputPrefix + "noMatch2.fastq").c_str(), ios::out));
	boost::shared_ptr<ofstream> multiMatch1(
			new ofstream((outputPrefix + "multiMatch.fastq").c_str(),
					ios::out));
	boost::shared_ptr<ofstream> multiMatch2(
			new ofstream((outputPrefix + "multiMatch.fastq").c_str(),
					ios::out));
	outputFiles["noMatch1"] = noMatch1;
	outputFiles["noMatch2"] = noMatch2;
	outputFiles["multiMatch1"] = multiMatch1;
	outputFiles["multiMatch2"] = multiMatch2;

	//initialize variables and print filter ids
	for (vector<string>::const_iterator j = hashSigs.begin();
			j != hashSigs.end(); ++j)
	{
		const vector<string> idsInFilter = (*filters[*j]).getFilterIds();
		for (vector<string>::const_iterator i = idsInFilter.begin();
				i != idsInFilter.end(); ++i)
		{
			boost::shared_ptr<ofstream> temp1(
					new ofstream((outputPrefix + *i + "2.fastq").c_str(),
							ios::out));
			boost::shared_ptr<ofstream> temp2(
					new ofstream((outputPrefix + *i + "2.fastq").c_str(),
							ios::out));
			outputFiles[*i + "1"] = temp1;
			outputFiles[*i + "2"] = temp2;
			readStatusOutput << "\t" << *i << "_"
					<< (*(infoFiles[*j].front())).getKmerSize();
			aboveThreshold[*i + "1"] = 0;
			belowThreshold[*i + "2"] = 0;
		}
	}
	readStatusOutput << "\n";

	//Todo: make sure this prints out only when filters are loaded
	//gcc currently optimizes to print this before loading can complete
	cerr << "Filtering Start" << "\n";

	FastaReader sequence1(file1.c_str(), FastaReader::NO_FOLD_CASE);
	FastaReader sequence2(file2.c_str(), FastaReader::NO_FOLD_CASE);
	FastqRecord rec1;
	FastqRecord rec2;
	//hits results stored in hashmap of filter names and hits
	boost::unordered_map<string, size_t> hits1(filterNum);
	boost::unordered_map<string, size_t> hits2(filterNum);

	while (sequence1 >> rec1 && sequence2 >> rec2) {
		//split reads into kmerSizes specified (ignore trailing bases)
		//for skipping bad reads
		bool readOK = true;

		//initialize hits to zero
		for (vector<string>::const_iterator j = hashSigs.begin();
				j != hashSigs.end(); ++j)
		{
			const vector<string> &idsInFilter = (*filters[*j]).getFilterIds();
			for (vector<string>::const_iterator i = idsInFilter.begin();
					i != idsInFilter.end(); ++i)
			{
				hits1[*i] = 0;
				hits2[*i] = 0;
			}
		}

		//for each hashSigniture/kmer combo multi, cut up read into kmer sized used
		for (vector<string>::const_iterator j = hashSigs.begin();
				j != hashSigs.end(); ++j)
		{
			string tempStr1 = rec1.id.substr(0, rec1.id.find_last_of("/"));
			string tempStr2 = rec2.id.substr(0, rec2.id.find_last_of("/"));
			if (tempStr1 == tempStr2) {
				if (!evaluateRead(rec1, *j, hits1)
						|| !evaluateRead(rec2, *j, hits2))
				{
					readOK = false;
					break;
				}
			} else {
				cerr << "Read IDs do not match" << "\n" << tempStr1 << "\n"
						<< tempStr2 << endl;
				exit(1);
			}
		}

		//print readID
		readStatusOutput << rec1.id << "\t" << rec1.seq.length();

		++totalReads;
		if (totalReads % 100000 == 0) {
			cout << "Currently Reading Read Number: " << totalReads << endl;
		}
		if (readOK) {
			int16_t totalHits = 0;
			for (vector<string>::const_iterator j = hashSigs.begin();
					j != hashSigs.end(); ++j)
			{
				//update summary
				const vector<string> &idsInFilter =
						(*filters[*j]).getFilterIds();
				for (vector<string>::const_iterator i = idsInFilter.begin();
						i != idsInFilter.end(); ++i)
				{
					//print read status
					readStatusOutput << "\t" << hits1[*i] << '/' << hits2[*i];

					//pick threshold, by percent or by absolute value
					int16_t kmerSize = (*(infoFiles[*j].front())).getKmerSize();
					size_t threshold1 = size_t(
							percentMinHit * (rec1.seq.length() / kmerSize));
					size_t threshold2 = size_t(
							percentMinHit * (rec2.seq.length() / kmerSize));
					if (minHit > threshold1) {
						threshold1 = minHit;
					}
					if (minHit > threshold2) {
						threshold2 = minHit;
					}

					if (hits1[*i] >= threshold1 && hits2[*i] >= threshold2) {
						++totalHits;
						++aboveThreshold[*i];
					} else if (hits1[*i] != 0 && hits2[*i] != 0) {
						++belowThreshold[*i];
					}
				}
			}
			if (totalHits == 0) {
				(*outputFiles["noMatch1"]) << "@" << rec1.id << "\n" << rec1.seq
						<< "\n+\n" << rec1.qual << "\n";
				(*outputFiles["noMatch2"]) << "@" << rec2.id << "\n" << rec2.seq
						<< "\n+\n" << rec2.qual << "\n";
			} else if (totalHits > 1) {
				(*outputFiles["multiMatch1"]) << "@" << rec1.id << "\n"
						<< rec1.seq << "\n+\n" << rec1.qual << "\n";
				(*outputFiles["multiMatch2"]) << "@" << rec2.id << "\n"
						<< rec2.seq << "\n+\n" << rec2.qual << "\n";
			} else {
				for (vector<string>::const_iterator j = hashSigs.begin();
						j != hashSigs.end(); ++j)
				{
					const vector<string> idsInFilter =
							(*filters[*j]).getFilterIds();
					for (vector<string>::const_iterator i = idsInFilter.begin();
							i != idsInFilter.end(); ++i)
					{
						//pick threshold, by percent or by absolute value
						int16_t kmerSize =
								(*(infoFiles[*j].front())).getKmerSize();
						size_t threshold1 = size_t(
								percentMinHit * (rec1.seq.length() / kmerSize));
						size_t threshold2 = size_t(
								percentMinHit * (rec2.seq.length() / kmerSize));
						if (minHit > threshold1) {
							threshold1 = minHit;
						}
						if (minHit > threshold2) {
							threshold2 = minHit;
						}

						if (hits1[*i] >= threshold1 && hits2[*i] >= threshold2)
						{
							(*outputFiles[*i + "1"]) << "@" << rec1.id << "\n"
									<< rec1.seq << "\n+\n" << rec1.qual << "\n";
							(*outputFiles[*i + "2"]) << "@" << rec2.id << "\n"
									<< rec2.seq << "\n+\n" << rec2.qual << "\n";
							break;
						}
					}
				}
			}
		} else {
			for (vector<string>::const_iterator j = hashSigs.begin();
					j != hashSigs.end(); ++j)
			{
				const vector<string> &idsInFilter =
						(*filters[*j]).getFilterIds();
				for (vector<string>::const_iterator i = idsInFilter.begin();
						i != idsInFilter.end(); ++i)
				{
					//print read status
					readStatusOutput << "\t" << "na/na";
				}
			}
			(*outputFiles["noMatch1"]) << "@" << rec1.id << "\n" << rec1.seq
					<< "\n+\n" << rec1.qual << "\n";
			(*outputFiles["noMatch2"]) << "@" << rec2.id << "\n" << rec2.seq
					<< "\n+\n" << rec2.qual << "\n";

		}
		readStatusOutput << "\n";
	}
	if(!(sequence2 >> rec2) && sequence1.eof() && sequence2.eof())
	{
		cerr << "error: eof bit not flipped. Input files may be different lengths" << endl;
	}

	//close sorting files
	for (vector<string>::const_iterator j = hashSigs.begin();
			j != hashSigs.end(); ++j)
	{
		const vector<string> idsInFilter = (*filters[*j]).getFilterIds();
		for (vector<string>::const_iterator i = idsInFilter.begin();
				i != idsInFilter.end(); ++i)
		{
			outputFiles[*i + "1"]->flush();
			outputFiles[*i + "2"]->close();
			outputFiles[*i + "1"]->flush();
			outputFiles[*i + "2"]->close();
		}
	}
	outputFiles["noMatch1"]->flush();
	outputFiles["noMatch1"]->close();
	outputFiles["noMatch2"]->flush();
	outputFiles["noMatch2"]->close();
	outputFiles["multiMatch1"]->flush();
	outputFiles["multiMatch1"]->close();
	outputFiles["multiMatch2"]->flush();
	outputFiles["multiMatch2"]->close();
	cout << "Total Reads:" << totalReads << endl;
	printSummary(outputPrefix, aboveThreshold, belowThreshold, totalReads);
}

/*
 * Filters reads -> uses paired end information
 * Assumes only one hash signature exists (load only filters with same
 * hash functions)
 */
void BioBloomClassifier::filterPair(const string &file1, const string &file2,
		const string &outputPrefix)
{
	//check if output folder exists
	string tempStr = outputPrefix.substr(0, outputPrefix.find_last_of("/"));
	folderCheck(tempStr);

	ofstream readStatusOutput((outputPrefix + "_status.tsv").c_str(), ios::out);

	//print header
	readStatusOutput << "readID\tseqSize";

	//variables for storing results summary
	boost::unordered_map<string, size_t> aboveThreshold;
	boost::unordered_map<string, size_t> belowThreshold;
	size_t totalReads = 0;

	//initialize variables and print filter ids
	for (vector<string>::const_iterator j = hashSigs.begin();
			j != hashSigs.end(); ++j)
	{
		const vector<string> idsInFilter = (*filters[*j]).getFilterIds();
		for (vector<string>::const_iterator i = idsInFilter.begin();
				i != idsInFilter.end(); ++i)
		{
			readStatusOutput << "\t" << *i << "_"
					<< (*(infoFiles[*j].front())).getKmerSize();
			aboveThreshold[*i + "1"] = 0;
			belowThreshold[*i + "2"] = 0;
		}
	}
	readStatusOutput << "\n";

	//Todo: make sure this prints out only when filters are loaded
	//gcc currently optimizes to print this before loading can complete
	cerr << "Filtering Start" << "\n";

	FastaReader sequence1(file1.c_str(), FastaReader::NO_FOLD_CASE);
	FastaReader sequence2(file2.c_str(), FastaReader::NO_FOLD_CASE);
	FastqRecord rec1;
	FastqRecord rec2;
	//hits results stored in hashmap of filter names and hits
	boost::unordered_map<string, size_t> hits1(filterNum);
	boost::unordered_map<string, size_t> hits2(filterNum);

	while (sequence1 >> rec1 && sequence2 >> rec2) {
		//split reads into kmerSizes specified (ignore trailing bases)
		//for skipping bad reads
		bool readOK = true;

		//initialize hits to zero
		for (vector<string>::const_iterator j = hashSigs.begin();
				j != hashSigs.end(); ++j)
		{
			const vector<string> &idsInFilter = (*filters[*j]).getFilterIds();
			for (vector<string>::const_iterator i = idsInFilter.begin();
					i != idsInFilter.end(); ++i)
			{
				hits1[*i] = 0;
				hits2[*i] = 0;
			}
		}

		//for each hashSigniture/kmer combo multi, cut up read into kmer sized used
		for (vector<string>::const_iterator j = hashSigs.begin();
				j != hashSigs.end(); ++j)
		{
			string tempStr1 = rec1.id.substr(0, rec1.id.find_last_of("/"));
			string tempStr2 = rec2.id.substr(0, rec2.id.find_last_of("/"));
			if (tempStr1 == tempStr2) {
				if (!evaluateRead(rec1, *j, hits1)
						|| !evaluateRead(rec2, *j, hits2))
				{
					readOK = false;
					break;
				}
			} else {
				cerr << "Read IDs do not match" << "\n" << tempStr1 << "\n"
						<< tempStr2 << endl;
				exit(1);
			}
		}

		//print readID
		readStatusOutput << rec1.id << "\t" << rec1.seq.length();

		++totalReads;
		if (totalReads % 100000 == 0) {
			cout << "Currently Reading Read Number: " << totalReads << endl;
		}
		if (readOK) {
			for (vector<string>::const_iterator j = hashSigs.begin();
					j != hashSigs.end(); ++j)
			{
				//update summary
				const vector<string> &idsInFilter =
						(*filters[*j]).getFilterIds();
				for (vector<string>::const_iterator i = idsInFilter.begin();
						i != idsInFilter.end(); ++i)
				{
					//print read status
					readStatusOutput << "\t" << hits1[*i] << '/' << hits2[*i];

					//pick threshold, by percent or by absolute value
					int16_t kmerSize = (*(infoFiles[*j].front())).getKmerSize();
					size_t threshold1 = size_t(
							percentMinHit * (rec1.seq.length() / kmerSize));
					size_t threshold2 = size_t(
							percentMinHit * (rec2.seq.length() / kmerSize));
					if (minHit > threshold1) {
						threshold1 = minHit;
					}
					if (minHit > threshold2) {
						threshold2 = minHit;
					}

					if (hits1[*i] >= threshold1 && hits2[*i] >= threshold2) {
						++aboveThreshold[*i];
					} else if (hits1[*i] != 0 && hits2[*i] != 0) {
						++belowThreshold[*i];
					}
				}
			}
		} else {
			for (vector<string>::const_iterator j = hashSigs.begin();
					j != hashSigs.end(); ++j)
			{
				const vector<string> &idsInFilter =
						(*filters[*j]).getFilterIds();
				for (vector<string>::const_iterator i = idsInFilter.begin();
						i != idsInFilter.end(); ++i)
				{
					//print read status
					readStatusOutput << "\t" << "na/na";
				}
			}
		}
		readStatusOutput << "\n";
	}
	if(!(sequence2 >> rec2) && sequence1.eof() && sequence2.eof())
	{
		cerr << "error: eof bit not flipped. Input files may be different lengths" << endl;
	}
	cout << "Total Reads:" << totalReads << endl;
	printSummary(outputPrefix, aboveThreshold, belowThreshold, totalReads);
}

/*
 * For a single read evaluate hits for a single hash signature
 * Returns true if read is valid (no ambiguity characters)
 * Updates hits value to number of hits (hashSig is used to as key)
 */
bool BioBloomClassifier::evaluateRead(const FastqRecord &rec,
		const string &hashSig, boost::unordered_map<string, size_t> &hits)
{
	bool readOK = true;

	//get filterIDs to iterate through has in a consistent order
	const vector<string> &idsInFilter = (*filters[hashSig]).getFilterIds();

	//get kmersize for set of info files
	int16_t kmerSize = (*(infoFiles[hashSig].front())).getKmerSize();

	//Establish tiling pattern
	int16_t startModifier1 = (rec.seq.length() % kmerSize) / 2;
	size_t currentKmerNum = 0;

	ReadsProcessor proc(kmerSize);
	//cut read into kmer size given
	while (rec.seq.length() >= (currentKmerNum + 1) * kmerSize) {

		const string &currentKmer = proc.prepSeq(rec.seq,
				currentKmerNum * kmerSize + startModifier1);

		//check to see if string is invalid
		if (!currentKmer.empty()) {
			const boost::unordered_map<string, bool> &results =
					filters[hashSig]->multiContains(currentKmer);

			//record hit number in order
			for (vector<string>::const_iterator i = idsInFilter.begin();
					i != idsInFilter.end(); ++i)
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
	return readOK;
}

//helper methods

/*
 * checks if file exists
 */
bool BioBloomClassifier::fexists(const string &filename) const
{
	ifstream ifile(filename.c_str());
	return ifile;
}

/*
 * checks if folder exists
 */
void BioBloomClassifier::folderCheck(const string &path) const
{
	struct stat sb;

	if (stat(path.c_str(), &sb) == 0) {
		if (!S_ISDIR(sb.st_mode)) {
			cerr << "Output folder: File exists with this name. " << path
					<< endl;
			exit(1);
		}
	} else {
		cerr << "Output folder: Does not exist. " << path << endl;
		exit(1);
	}
}

BioBloomClassifier::~BioBloomClassifier()
{
}

