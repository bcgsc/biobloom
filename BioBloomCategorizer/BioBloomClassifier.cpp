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
#include "ResultsManager.h"
#include "Common/Options.h"
#include <map>
#if _OPENMP
# include <omp.h>
#endif

BioBloomClassifier::BioBloomClassifier(const vector<string> &filterFilePaths,
		double scoreThreshold, const string &prefix,
		const string &outputPostFix, unsigned minHit, bool minHitOnly) :
		m_scoreThreshold(scoreThreshold), m_filterNum(filterFilePaths.size()), m_prefix(
				prefix), m_postfix(outputPostFix), m_minHit(minHit), m_mode(
				STD), m_mainFilter(""), m_inclusive(false)
{
	if (minHitOnly) {
		m_mode = MINHITONLY;
	}
	loadFilters(filterFilePaths);
	if (m_scoreThreshold == 1) {
		m_mode = BESTHIT;
		//TODO: make it possible later?
		//best hit will not allow for more than one hash sig.
		assert(m_mode == BESTHIT && m_hashSigs.size() == 1);
	}
}

/*
 * Generic filtering function (single end, no fa or fq file outputs)
 */
void BioBloomClassifier::filter(const vector<string> &inputFiles)
{

	//results summary object
	ResultsManager resSummary(m_hashSigs, m_filters, m_infoFiles, m_inclusive);

	size_t totalReads = 0;

	//print out header info and initialize variables

	cerr << "Filtering Start" << endl;

	for (vector<string>::const_iterator it = inputFiles.begin();
			it != inputFiles.end(); ++it)
	{
		FastaReader sequence(it->c_str(), FastaReader::NO_FOLD_CASE);
#pragma omp parallel
		for (FastqRecord rec;;) {
			bool good;
#pragma omp critical(sequence)
			{
				good = sequence >> rec;
				//track read progress
			}
			if (good) {
#pragma omp critical(totalReads)
				{
					++totalReads;
					if (totalReads % 10000000 == 0) {
						cerr << "Currently Reading Read Number: " << totalReads
								<< endl;
					}
				}
				unordered_map<string, bool> hits(m_filterNum);
				double score = 0; //Todo: figure out what happens to this if multiple hashSigs are used

				//for each hashSigniture/kmer combo multi, cut up read into kmer sized used
				for (vector<string>::const_iterator j = m_hashSigs.begin();
						j != m_hashSigs.end(); ++j)
				{
					evaluateRead(rec, *j, hits, score);
				}

				//Evaluate hit data and record for summary and print if needed
				printSingle(rec, score, resSummary.updateSummaryData(hits));

			} else
				break;
		}
		assert(sequence.eof());
	}

	cerr << "Total Reads:" << totalReads << endl;

	cerr << "Writing file: " << m_prefix + "_summary.tsv" << endl;

	Dynamicofstream summaryOutput(m_prefix + "_summary.tsv");
	summaryOutput << resSummary.getResultsSummary(totalReads);
	summaryOutput.close();
}

/*
 * Filters reads
 * Assumes only one hash signature exists (load only filters with same
 * hash functions)
 * Prints reads into seperate files
 */
void BioBloomClassifier::filterPrint(const vector<string> &inputFiles,
		const string &outputType)
{

	//results summary object
	ResultsManager resSummary(m_hashSigs, m_filters, m_infoFiles, m_inclusive);

	size_t totalReads = 0;

	unordered_map<string, boost::shared_ptr<Dynamicofstream> > outputFiles;
	boost::shared_ptr<Dynamicofstream> no_match(
			new Dynamicofstream(
					m_prefix + "_" + NO_MATCH + "." + outputType + m_postfix));
	boost::shared_ptr<Dynamicofstream> multi_match(
			new Dynamicofstream(
					m_prefix + "_" + MULTI_MATCH + "." + outputType
							+ m_postfix));
	outputFiles[NO_MATCH] = no_match;
	outputFiles[MULTI_MATCH] = multi_match;

	//initialize variables
	for (vector<string>::const_iterator j = m_hashSigs.begin();
			j != m_hashSigs.end(); ++j)
	{
		const vector<string> idsInFilter = (*m_filters[*j]).getFilterIds();
		for (vector<string>::const_iterator i = idsInFilter.begin();
				i != idsInFilter.end(); ++i)
		{
			boost::shared_ptr<Dynamicofstream> temp(
					new Dynamicofstream(
							m_prefix + "_" + *i + "." + outputType
									+ m_postfix));
			outputFiles[*i] = temp;
		}
	}

	//print out header info and initialize variables

	cerr << "Filtering Start" << endl;

	for (vector<string>::const_iterator it = inputFiles.begin();
			it != inputFiles.end(); ++it)
	{
		FastaReader sequence(it->c_str(), FastaReader::NO_FOLD_CASE);
#pragma omp parallel
		for (FastqRecord rec;;) {
			bool good;
#pragma omp critical(sequence)
			{
				good = sequence >> rec;
			}
			if (good) {
#pragma omp critical(totalReads)
				{
					++totalReads;
					if (totalReads % 10000000 == 0) {
						cerr << "Currently Reading Read Number: " << totalReads
								<< endl;
					}
				}
				unordered_map<string, bool> hits(m_filterNum);
				double score = 0.0; //Todo: figure out what happens to this if multiple hashSigs are used

				//for each hashSigniture/kmer combo multi, cut up read into kmer sized used
				for (vector<string>::const_iterator j = m_hashSigs.begin();
						j != m_hashSigs.end(); ++j)
				{
					evaluateRead(rec, *j, hits, score);
				}

				//Evaluate hit data and record for summary
				const string &outputFileName = resSummary.updateSummaryData(
						hits);

				printSingle(rec, score, outputFileName);

				printSingleToFile(outputFileName, rec, outputFiles, outputType,
						score);

			} else
				break;
		}
		assert(sequence.eof());
	}

	//close sorting files
	for (unordered_map<string, boost::shared_ptr<Dynamicofstream> >::iterator j =
			outputFiles.begin(); j != outputFiles.end(); ++j)
	{
		j->second->close();
		cerr << "File written to: "
				<< m_prefix + "_" + j->first + "." + outputType + m_postfix
				<< endl;
	}
	cerr << "Total Reads:" << totalReads << endl;
	cerr << "Writing file: " << m_prefix + "_summary.tsv" << endl;

	Dynamicofstream summaryOutput(m_prefix + "_summary.tsv");
	summaryOutput << resSummary.getResultsSummary(totalReads);
	summaryOutput.close();
}

/*
 * Filters reads -> uses paired end information
 * Assumes only one hash signature exists (load only filters with same
 * hash functions)
 */
void BioBloomClassifier::filterPair(const string &file1, const string &file2)
{

	//results summary object
	ResultsManager resSummary(m_hashSigs, m_filters, m_infoFiles, m_inclusive);

	size_t totalReads = 0;

	cerr << "Filtering Start" << "\n";

	FastaReader sequence1(file1.c_str(), FastaReader::NO_FOLD_CASE);
	FastaReader sequence2(file2.c_str(), FastaReader::NO_FOLD_CASE);
#pragma omp parallel
	for (FastqRecord rec1;;) {
		FastqRecord rec2;
		bool good1;
		bool good2;

#pragma omp critical(sequence1)
		{
			good1 = sequence1 >> rec1;
			good2 = sequence2 >> rec2;
			//track read progress
		}

		if (good1 && good2) {
#pragma omp critical(totalReads)
			{
				++totalReads;
				if (totalReads % 10000000 == 0) {
					cerr << "Currently Reading Read Number: " << totalReads
							<< endl;
				}
			}

			//hits results stored in hashmap of filter names and hits
			unordered_map<string, bool> hits1(m_filterNum);
			unordered_map<string, bool> hits2(m_filterNum);

			double score1;
			double score2;

			//for each hashSigniture/kmer combo multi, cut up read into kmer sized used
			for (vector<string>::const_iterator j = m_hashSigs.begin();
					j != m_hashSigs.end(); ++j)
			{
				string tempStr1 = rec1.id.substr(0, rec1.id.find_last_of("/"));
				string tempStr2 = rec2.id.substr(0, rec2.id.find_last_of("/"));
				if (tempStr1 == tempStr2) {
					evaluateRead(rec1, *j, hits1, score1);
					evaluateRead(rec2, *j, hits2, score2);
				} else {
					cerr << "Read IDs do not match" << "\n" << tempStr1 << "\n"
							<< tempStr2 << endl;
					exit(1);
				}
			}

			//Evaluate hit data and record for summary and print if needed
			printPair(rec1, rec2, score1, score2,
					resSummary.updateSummaryData(hits1, hits2));
		} else
			break;
	}
	if (!sequence1.eof() || !sequence2.eof()) {
		cerr
				<< "error: eof bit not flipped. Input files may be different lengths"
				<< endl;
	}

	cerr << "Total Reads:" << totalReads << endl;
	cerr << "Writing file: " << m_prefix + "_summary.tsv" << endl;

	Dynamicofstream summaryOutput(m_prefix + "_summary.tsv");
	summaryOutput << resSummary.getResultsSummary(totalReads);
	summaryOutput.close();
}

/*
 * Filters reads -> uses paired end information
 * Assumes only one hash signature exists (load only filters with same
 * hash functions)
 * prints reads
 */
void BioBloomClassifier::filterPairPrint(const string &file1,
		const string &file2, const string &outputType)
{

	//results summary object
	ResultsManager resSummary(m_hashSigs, m_filters, m_infoFiles, m_inclusive);

	size_t totalReads = 0;

	unordered_map<string, boost::shared_ptr<Dynamicofstream> > outputFiles;
	boost::shared_ptr<Dynamicofstream> noMatch1(
			new Dynamicofstream(
					m_prefix + "_" + NO_MATCH + "_1." + outputType
							+ m_postfix));
	boost::shared_ptr<Dynamicofstream> noMatch2(
			new Dynamicofstream(
					m_prefix + "_" + NO_MATCH + "_2." + outputType
							+ m_postfix));
	boost::shared_ptr<Dynamicofstream> multiMatch1(
			new Dynamicofstream(
					m_prefix + "_" + MULTI_MATCH + "_1." + outputType
							+ m_postfix));
	boost::shared_ptr<Dynamicofstream> multiMatch2(
			new Dynamicofstream(
					m_prefix + "_" + MULTI_MATCH + "_2." + outputType
							+ m_postfix));
	outputFiles[NO_MATCH + "1"] = noMatch1;
	outputFiles[NO_MATCH + "2"] = noMatch2;
	outputFiles[MULTI_MATCH + "1"] = multiMatch1;
	outputFiles[MULTI_MATCH + "2"] = multiMatch2;

	//initialize variables and print filter ids
	for (vector<string>::const_iterator j = m_hashSigs.begin();
			j != m_hashSigs.end(); ++j)
	{
		const vector<string> idsInFilter = (*m_filters[*j]).getFilterIds();
		for (vector<string>::const_iterator i = idsInFilter.begin();
				i != idsInFilter.end(); ++i)
		{
			boost::shared_ptr<Dynamicofstream> temp1(
					new Dynamicofstream(
							m_prefix + "_" + *i + "_1." + outputType
									+ m_postfix));
			boost::shared_ptr<Dynamicofstream> temp2(
					new Dynamicofstream(
							m_prefix + "_" + *i + "_2." + outputType
									+ m_postfix));
			outputFiles[*i + "1"] = temp1;
			outputFiles[*i + "2"] = temp2;
		}
	}

	//for output files in consistent order

	cerr << "Filtering Start" << "\n";

	FastaReader sequence1(file1.c_str(), FastaReader::NO_FOLD_CASE);
	FastaReader sequence2(file2.c_str(), FastaReader::NO_FOLD_CASE);
#pragma omp parallel
	for (FastqRecord rec1;;) {
		FastqRecord rec2;
		bool good1;
		bool good2;

#pragma omp critical(sequence1)
		{
			good1 = sequence1 >> rec1;
			good2 = sequence2 >> rec2;
		}

		if (good1 && good2) {
#pragma omp critical(totalReads)
			{
				++totalReads;
				if (totalReads % 10000000 == 0) {
					cerr << "Currently Reading Read Number: " << totalReads
							<< endl;
				}
			}

			//hits results stored in hashmap of filter names and hits
			unordered_map<string, bool> hits1(m_filterNum);
			unordered_map<string, bool> hits2(m_filterNum);

			double score1;
			double score2;
			//for each hashSigniture/kmer combo multi, cut up read into kmer sized used
			for (vector<string>::const_iterator j = m_hashSigs.begin();
					j != m_hashSigs.end(); ++j)
			{
				string tempStr1 = rec1.id.substr(0, rec1.id.find_last_of("/"));
				string tempStr2 = rec2.id.substr(0, rec2.id.find_last_of("/"));
				if (tempStr1 == tempStr2) {
					evaluateRead(rec1, *j, hits1, score1);
					evaluateRead(rec2, *j, hits2, score2);
				} else {
					cerr << "Read IDs do not match" << "\n" << tempStr1 << "\n"
							<< tempStr2 << endl;
					exit(1);
				}
			}

			//Evaluate hit data and record for summary

			const string &outputFileName = resSummary.updateSummaryData(hits1,
					hits2);
			printPair(rec1, rec2, score1, score2, outputFileName);
			printPairToFile(outputFileName, rec1, rec2, outputFiles, outputType,
					score1, score2);
		} else
			break;
	}
	if (!sequence1.eof() || !sequence2.eof()) {
		cerr
				<< "error: eof bit not flipped. Input files may be different lengths"
				<< endl;
	}

	//close sorting files
	for (unordered_map<string, boost::shared_ptr<Dynamicofstream> >::iterator j =
			outputFiles.begin(); j != outputFiles.end(); ++j)
	{
		j->second->close();
		cerr << "File written to: "
				<< m_prefix + "_" + j->first + "." + outputType + m_postfix
				<< endl;
	}

	cerr << "Total Reads:" << totalReads << endl;
	cerr << "Writing file: " << m_prefix + "_summary.tsv" << endl;

	Dynamicofstream summaryOutput(m_prefix + "_summary.tsv");
	summaryOutput << resSummary.getResultsSummary(totalReads);
	summaryOutput.close();
}

/*
 * Filters reads -> uses paired end information
 * Assumes only one hash signature exists (load only filters with same
 * hash functions)
 */
void BioBloomClassifier::filterPairBAM(const string &file)
{

	//results summary object
	ResultsManager resSummary(m_hashSigs, m_filters, m_infoFiles, m_inclusive);

	unordered_map<string, FastqRecord> unPairedReads;

	size_t totalReads = 0;

	//print out header info and initialize variables for summary

	cerr << "Filtering Start" << "\n";

	FastaReader sequence(file.c_str(), FastaReader::NO_FOLD_CASE);
#pragma omp parallel
	for (FastqRecord rec;;) {
		bool good;
		bool pairfound;
		FastqRecord rec1;
		FastqRecord rec2;
#pragma omp critical(unPairedReads)
		{
			good = sequence >> rec;
			//track read progress
			if (good) {
				string readID = rec.id.substr(0, rec.id.length() - 2);
				if (unPairedReads.find(readID) != unPairedReads.end()) {
					rec1 = rec.id.at(rec.id.length() - 1) == '1' ?
							rec : unPairedReads[readID];
					rec2 = rec.id.at(rec.id.length() - 1) == '2' ?
							rec : unPairedReads[readID];
					pairfound = true;
					unPairedReads.erase(readID);
				} else {
					pairfound = false;
					unPairedReads[readID] = rec;
				}
			}
		}

		if (good) {
			if (pairfound) {
#pragma omp critical(totalReads)
				{
					++totalReads;
					if (totalReads % 10000000 == 0) {
						cerr << "Currently Reading Read Number: " << totalReads
								<< endl;
					}
				}

				unordered_map<string, bool> hits1(m_filterNum);
				unordered_map<string, bool> hits2(m_filterNum);

				double score1;
				double score2;
				//for each hashSigniture/kmer combo multi, cut up read into kmer sized used
				for (vector<string>::const_iterator j = m_hashSigs.begin();
						j != m_hashSigs.end(); ++j)
				{
					evaluateRead(rec1, *j, hits1, score1);
					evaluateRead(rec2, *j, hits2, score2);
				}

				//Evaluate hit data and record for summary
				printPair(rec1, rec2, score1, score2,
						resSummary.updateSummaryData(hits1, hits2));
			}
		} else
			break;
	}
	assert(sequence.eof());

	cerr << "Total Reads:" << totalReads << endl;
	cerr << "Writing file: " << m_prefix + "_summary.tsv" << endl;

	Dynamicofstream summaryOutput(m_prefix + "_summary.tsv");
	summaryOutput << resSummary.getResultsSummary(totalReads);
	summaryOutput.close();
}

/*
 * Filters reads -> uses paired end information
 * Assumes only one hash signature exists (load only filters with same
 * hash functions)
 * Prints reads into separate files
 */
void BioBloomClassifier::filterPairBAMPrint(const string &file,
		const string &outputType)
{

	//results summary object
	ResultsManager resSummary(m_hashSigs, m_filters, m_infoFiles, m_inclusive);

	unordered_map<string, FastqRecord> unPairedReads;

	size_t totalReads = 0;

	unordered_map<string, boost::shared_ptr<Dynamicofstream> > outputFiles;
	boost::shared_ptr<Dynamicofstream> noMatch1(
			new Dynamicofstream(
					m_prefix + "_" + NO_MATCH + "_1." + outputType
							+ m_postfix));
	boost::shared_ptr<Dynamicofstream> noMatch2(
			new Dynamicofstream(
					m_prefix + "_" + NO_MATCH + "_2." + outputType
							+ m_postfix));
	boost::shared_ptr<Dynamicofstream> multiMatch1(
			new Dynamicofstream(
					m_prefix + "_" + MULTI_MATCH + "_1." + outputType
							+ m_postfix));
	boost::shared_ptr<Dynamicofstream> multiMatch2(
			new Dynamicofstream(
					m_prefix + "_" + MULTI_MATCH + "_2." + outputType
							+ m_postfix));
	outputFiles[NO_MATCH + "1"] = noMatch1;
	outputFiles[NO_MATCH + "2"] = noMatch2;
	outputFiles[MULTI_MATCH + "1"] = multiMatch1;
	outputFiles[MULTI_MATCH + "2"] = multiMatch2;

	//initialize variables and print filter ids
	for (vector<string>::const_iterator j = m_hashSigs.begin();
			j != m_hashSigs.end(); ++j)
	{
		const vector<string> idsInFilter = (*m_filters[*j]).getFilterIds();
		for (vector<string>::const_iterator i = idsInFilter.begin();
				i != idsInFilter.end(); ++i)
		{
			boost::shared_ptr<Dynamicofstream> temp1(
					new Dynamicofstream(
							m_prefix + "_" + *i + "_1." + outputType
									+ m_postfix));
			boost::shared_ptr<Dynamicofstream> temp2(
					new Dynamicofstream(
							m_prefix + "_" + *i + "_2." + outputType
									+ m_postfix));
			outputFiles[*i + "1"] = temp1;
			outputFiles[*i + "2"] = temp2;
		}
	}

	//print out header info and initialize variables for summary
	cerr << "Filtering Start" << "\n";

	FastaReader sequence(file.c_str(), FastaReader::NO_FOLD_CASE);
#pragma omp parallel
	for (FastqRecord rec;;) {
		bool good;
		bool pairfound = false;
		FastqRecord rec1;
		FastqRecord rec2;
#pragma omp critical(unPairedReads)
		{
			good = sequence >> rec;
			//track read progress
			if (good) {
				string readID = rec.id.substr(0, rec.id.length() - 2);
				if (unPairedReads.find(readID) != unPairedReads.end()) {
					rec1 = rec.id.at(rec.id.length() - 1) == '1' ?
							rec : unPairedReads[readID];
					rec2 = rec.id.at(rec.id.length() - 1) == '2' ?
							rec : unPairedReads[readID];
					pairfound = true;
					unPairedReads.erase(readID);
				} else {
					unPairedReads[readID] = rec;
				}
			}
		}
		if (good) {
			if (pairfound) {
#pragma omp critical(totalReads)
				{
					++totalReads;
					if (totalReads % 10000000 == 0) {
						cerr << "Currently Reading Read Number: " << totalReads
								<< endl;
					}
				}

				unordered_map<string, bool> hits1(m_filterNum);
				unordered_map<string, bool> hits2(m_filterNum);

				double score1;
				double score2;

				//for each hashSigniture/kmer combo multi, cut up read into kmer sized used
				for (vector<string>::const_iterator j = m_hashSigs.begin();
						j != m_hashSigs.end(); ++j)
				{
					string tempStr1 = rec1.id.substr(0,
							rec1.id.find_last_of("/"));
					string tempStr2 = rec2.id.substr(0,
							rec2.id.find_last_of("/"));
					if (tempStr1 == tempStr2) {
						evaluateRead(rec1, *j, hits1, score1);
						evaluateRead(rec2, *j, hits2, score2);
					} else {
						cerr << "Read IDs do not match" << "\n" << tempStr1
								<< "\n" << tempStr2 << endl;
						exit(1);
					}
				}

				//Evaluate hit data and record for summary
				const string &outputFileName = resSummary.updateSummaryData(
						hits1, hits2);
				printPairToFile(outputFileName, rec1, rec2, outputFiles,
						outputType, score1, score2);
				printPair(rec1, rec2, score1, score2, outputFileName);
			}
		} else
			break;
	}
	assert(sequence.eof());

	//close sorting files
	for (unordered_map<string, boost::shared_ptr<Dynamicofstream> >::iterator j =
			outputFiles.begin(); j != outputFiles.end(); ++j)
	{
		j->second->close();
		cerr << "File written to: "
				<< m_prefix + "_" + j->first + "." + outputType + m_postfix
				<< endl;
	}

	cerr << "Total Reads:" << totalReads << endl;
	cerr << "Writing file: " << m_prefix + "_summary.tsv" << endl;

	Dynamicofstream summaryOutput(m_prefix + "_summary.tsv");
	summaryOutput << resSummary.getResultsSummary(totalReads);
	summaryOutput.close();
}

//helper methods

/*
 * Loads list of filters into memory
 * todo: Implement non-block I/O when loading multiple filters at once
 */
void BioBloomClassifier::loadFilters(const vector<string> &filterFilePaths)
{
	cerr << "Starting to Load Filters." << endl;
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
		hashSig << info->getHashNum() << info->getKmerSize();

		//if hashSig exists add filter to list
		if (m_infoFiles.count(hashSig.str()) != 1) {
			m_hashSigs.push_back(hashSig.str());
			vector<boost::shared_ptr<BloomFilterInfo> > tempVect;
			boost::shared_ptr<MultiFilter> temp(
					new MultiFilter(info->getHashNum(), info->getKmerSize()));
			m_filters[hashSig.str()] = temp;
			m_infoFiles[hashSig.str()] = tempVect;
		}
		m_infoFiles[hashSig.str()].push_back(info);
		boost::shared_ptr<BloomFilter> filter(
				new BloomFilter(info->getCalcuatedFilterSize(),
						info->getHashNum(), info->getKmerSize(), *it));
		m_filters[hashSig.str()]->addFilter(info->getFilterID(), filter);
		m_filtersSingle[info->getFilterID()] = filter;
		m_filterOrder.push_back(info->getFilterID());
		cerr << "Loaded Filter: " + info->getFilterID() << endl;
	}
	if (m_scoreThreshold == 1 && m_hashSigs.size() > 1) {
		cerr
				<< "If -s = 1 (best hit mode) all filters must use the same k and same number of hash functions."
				<< endl;
		exit(1);
	}
	cerr << "Filter Loading Complete." << endl;
}

/*
 * checks if file exists
 */
bool BioBloomClassifier::fexists(const string &filename) const
{
	ifstream ifile(filename.c_str());
	return ifile;
}

/*
 * Collaborative filtering method
 * Assume filters use the same k-mer size
 */
void BioBloomClassifier::evaluateReadCollab(const FastqRecord &rec,
		const string &hashSig, unordered_map<string, bool> &hits)
{
	//get filterIDs to iterate through has in a consistent order
	unsigned kmerSize = m_infoFiles.at(hashSig).front()->getKmerSize();

	//todo: read proc possibly unneeded, see evalSingle
	ReadsProcessor proc(kmerSize);

	//create storage for hits per filter
	std::multimap<unsigned, string> firstPassHits;

	//base for each filter until one filter obtains hit threshold
	//TODO: staggered pattering
	for (vector<string>::const_iterator i = m_filterOrder.begin();
			i != m_filterOrder.end(); ++i)
	{
		hits[*i] = false;
		unsigned screeningHits = 0;
		size_t screeningLoc = rec.seq.length() % kmerSize / 2;
		//First pass filtering
		while (rec.seq.length() >= screeningLoc + kmerSize) {
			const unsigned char* currentKmer = proc.prepSeq(rec.seq,
					screeningLoc);
			if (currentKmer != NULL) {
				if (m_filtersSingle.at(*i)->contains(currentKmer)) {
					++screeningHits;
				}
			}
			screeningLoc += kmerSize;
		}
		firstPassHits.insert(pair<unsigned, string>(screeningHits, *i));
	}

	double normalizationValue = rec.seq.length() - kmerSize + 1;
	double threshold = m_scoreThreshold * normalizationValue;
	size_t antiThreshold = static_cast<size_t>((1.0 - m_scoreThreshold) * normalizationValue);

	//evaluate promising group first
	for (multimap<unsigned, string>::reverse_iterator i =
			firstPassHits.rbegin(); i != firstPassHits.rend(); ++i)
	{
		string filterID = i->second;
		BloomFilter &tempFilter = *m_filtersSingle.at(filterID);
		if(SeqEval::evalSingle(rec, kmerSize, tempFilter, threshold, antiThreshold))
		{
			hits[filterID] = true;
			break;
		}
	}
}

/*
 * For a single read evaluate hits for a single hash signature
 * Sections with ambiguity bases are treated as misses
 * Updates hits value to number of hits (hashSig is used to as key)
 * Faster variant that assume there a redundant tile of 0
 */
void BioBloomClassifier::evaluateReadMin(const FastqRecord &rec,
		const string &hashSig, unordered_map<string, bool> &hits)
{
	//get filterIDs to iterate through has in a consistent order
	const vector<string> &idsInFilter = (*m_filters[hashSig]).getFilterIds();

	//get kmersize for set of info files
	unsigned kmerSize = m_infoFiles.at(hashSig).front()->getKmerSize();

	unordered_map<string, unsigned> tempHits;

	//Establish tiling pattern
	unsigned startModifier1 = (rec.seq.length() % kmerSize) / 2;
	size_t currentKmerNum = 0;

	for (vector<string>::const_iterator i = idsInFilter.begin();
			i != idsInFilter.end(); ++i)
	{
		tempHits[*i] = 0;
	}

	ReadsProcessor proc(kmerSize);
	//cut read into kmer size given
	while (rec.seq.length() >= (currentKmerNum + 1) * kmerSize) {

		const unsigned char* currentKmer = proc.prepSeq(rec.seq,
				currentKmerNum * kmerSize + startModifier1);

		//check to see if string is invalid
		if (currentKmer != NULL) {

			const unordered_map<string, bool> &results =
					m_filters[hashSig]->multiContains(currentKmer);

			//record hit number in order
			for (vector<string>::const_iterator i = idsInFilter.begin();
					i != idsInFilter.end(); ++i)
			{
				if (results.at(*i)) {
					++tempHits[*i];
				}
			}
		}
		++currentKmerNum;
	}
	for (vector<string>::const_iterator i = idsInFilter.begin();
			i != idsInFilter.end(); ++i)
	{
		hits[*i] = tempHits.at(*i) >= m_minHit;
	}
}

/*
 * For a single read evaluate hits for a single hash signature
 * Sections with ambiguity bases are treated as misses
 * Updates hits value to number of hits (hashSig is used to as key)
 */
void BioBloomClassifier::evaluateReadStd(const FastqRecord &rec,
		const string &hashSig, unordered_map<string, bool> &hits)
{

	//get filterIDs to iterate through has in a consistent order
	const vector<string> &idsInFilter = (*m_filters[hashSig]).getFilterIds();

	unsigned kmerSize = m_infoFiles.at(hashSig).front()->getKmerSize();

	//todo: read proc possibly unneeded, see evalSingle
	ReadsProcessor proc(kmerSize);

	double normalizationValue = rec.seq.length() - kmerSize + 1;
	double threshold = m_scoreThreshold * normalizationValue;
	size_t antiThreshold = static_cast<size_t>((1.0 - m_scoreThreshold) * normalizationValue);

	for (vector<string>::const_iterator i = idsInFilter.begin();
			i != idsInFilter.end(); ++i)
	{
		bool pass = false;
		hits[*i] = false;
		if (m_minHit > 0) {
			unsigned screeningHits = 0;
			size_t screeningLoc = rec.seq.length() % kmerSize / 2;
			//First pass filtering
			while (rec.seq.length() >= screeningLoc + kmerSize) {
				const unsigned char* currentKmer = proc.prepSeq(rec.seq,
						screeningLoc);
				if (currentKmer != NULL) {
					if (m_filtersSingle.at(*i)->contains(currentKmer)) {
						screeningHits++;
						if (screeningHits >= m_minHit) {
							pass = true;
							break;
						}
					}
				}
				screeningLoc += kmerSize;
			}
		} else {
			pass = true;
		}
		if (pass) {
			BloomFilter &tempFilter = *m_filtersSingle.at(*i);
			hits[*i] = SeqEval::evalSingle(rec, kmerSize, tempFilter, threshold, antiThreshold);
		}
	}
}

/*
 * For a single read evaluate hits for a single hash signature
 * Sections with ambiguity bases are treated as misses
 * Updates hits value to number of hits (hashSig is used to as key)
 */
double BioBloomClassifier::evaluateReadBestHit(const FastqRecord &rec,
		const string &hashSig, unordered_map<string, bool> &hits)
{
	//get filterIDs to iterate through has in a consistent order
	const vector<string> &idsInFilter = (*m_filters[hashSig]).getFilterIds();

	vector<string> bestFilters;
	double maxScore = 0;

	unsigned kmerSize = m_infoFiles.at(hashSig).front()->getKmerSize();

	ReadsProcessor proc(kmerSize);

	for (unsigned i = 0; i < idsInFilter.size(); ++i) {
		bool pass = false;
		hits[idsInFilter[i]] = false;
		if (m_minHit > 0) {
			unsigned screeningHits = 0;
			size_t screeningLoc = rec.seq.length() % kmerSize / 2;
			//First pass filtering
			while (rec.seq.length() >= screeningLoc + kmerSize) {
				const unsigned char* currentKmer = proc.prepSeq(rec.seq,
						screeningLoc);
				if (currentKmer != NULL) {
					if (m_filtersSingle.at(idsInFilter[i])->contains(
							currentKmer))
					{
						screeningHits++;
						if (screeningHits >= m_minHit) {
							pass = true;
							break;
						}
					}
				}
				screeningLoc += kmerSize;
			}
		} else {
			pass = true;
		}
		if (pass) {
			size_t currentLoc = 0;
			double score = 0;
			unsigned streak = 0;
			while (rec.seq.length() >= currentLoc + kmerSize) {
				const unsigned char* currentKmer = proc.prepSeq(rec.seq,
						currentLoc);
				if (streak == 0) {
					if (currentKmer != NULL) {
						if (m_filtersSingle.at(idsInFilter[i])->contains(
								currentKmer))
						{
							score += 0.5;
							++streak;
						}
						++currentLoc;
					} else {
						currentLoc += kmerSize + 1;
					}
				} else {
					if (currentKmer != NULL) {
						if (m_filtersSingle.at(idsInFilter[i])->contains(
								currentKmer))
						{
							++streak;
							score += 1 - 1 / (2 * streak);
							++currentLoc;
							continue;
						}
					} else {
						currentLoc += kmerSize + 1;
					}
					if (streak < opt::streakThreshold) {
						++currentLoc;
					} else {
						currentLoc += kmerSize;
					}
					streak = 0;
				}
			}
			if (maxScore < score) {
				maxScore = score;
				bestFilters.clear();
				bestFilters.push_back(idsInFilter[i]);
			} else if (maxScore == score) {
				bestFilters.push_back(idsInFilter[i]);
			}
		}
	}
	if (maxScore > 0) {
		for (unsigned i = 0; i < bestFilters.size(); ++i) {
			hits[bestFilters[i]] = true;
		}
	}
	return maxScore / (rec.seq.length() - kmerSize + 1);
}

void BioBloomClassifier::setMainFilter(const string &filtername)
{
	if (m_filtersSingle.find(filtername) == m_filtersSingle.end()) {
		cerr << "Filter with this name \"" << filtername
				<< "\" does not exist\n";
		cerr << "Valid filter Names:\n";
		for (vector<string>::const_iterator itr = m_filterOrder.begin();
				itr != m_filterOrder.end(); ++itr)
		{
			cerr << *itr << endl;
		}
		exit(1);
	}
	m_mainFilter = filtername;
}

BioBloomClassifier::~BioBloomClassifier()
{
}

