/*
 * BloomMapClassifier.cpp
 *
 *  Created on: Mar 23, 2016
 *      Author: cjustin
 */

#include <BloomMapClassifier.h>
#include <ResultsManager.hpp>
#include "Common/Options.h"
#include "Options.h"
#include <iostream>
#include <zlib.h>
#include "Common/kseq.h"
#include <boost/math/distributions/chi_squared.hpp>
//KSEQ_INIT(gzFile, gzread)

BloomMapClassifier::BloomMapClassifier(const string &filterFile) :
		m_filter(MIBloomFilter<ID>(filterFile)), m_numRead(0) {
	cerr << "FPR given allowed misses: " << m_filter.getFPR(opt::allowMisses)
			<< endl;
	//load in ID file
	string idFile = (filterFile).substr(0, (filterFile).length() - 3)
			+ "_ids.txt";
	google::dense_hash_map<ID, boost::shared_ptr<vector<ID> > > colliIDs;
	google::dense_hash_map<ID, string> ids;
	colliIDs.set_empty_key(opt::EMPTY);
	ids.set_empty_key(opt::EMPTY);
	if (!fexists(idFile)) {
		cerr << "Error: " << idFile.c_str()
				<< " File cannot be opened. A corresponding id file is needed."
				<< endl;
		exit(1);
	} else {
		cerr << "Loading ID file: " << idFile.c_str() << endl;
	}
	ifstream idFH(idFile.c_str(), ifstream::in);
	string line;
	getline(idFH, line);
	while (idFH.good()) {
		ID id;
		stringstream converter(line);
		converter >> id;

		colliIDs[id] = boost::shared_ptr<vector<ID> >(new vector<ID>());
		string fullID;
		unsigned idCount = 0;
		for (; converter >> fullID; ++idCount) {
			ID tempID = atoi(fullID.c_str());
			colliIDs[id]->push_back(tempID);
		}
		if (idCount == 1) {
			colliIDs[id]->clear();
			colliIDs[id]->push_back(id);
			ids[id] = fullID;
		}
		getline(idFH, line);
	}
	m_fullIDs = vector<string>(ids.size() + 1); //first element will be empty
	m_fullIDs[0] = "unknown"; //set first element to unknown
	m_idToIndex.set_empty_key("");
	for (google::dense_hash_map<ID, string>::const_iterator itr = ids.begin();
			itr != ids.end(); ++itr) {
//		cerr << itr->first << itr->second << endl;
		m_fullIDs[itr->first] = itr->second;
		m_idToIndex[itr->second] = itr->first;
	}
	m_colliIDs = vector<boost::shared_ptr<vector<ID> > >(colliIDs.size() + 1);
	for (google::dense_hash_map<ID, boost::shared_ptr<vector<ID> > >::const_iterator itr =
			colliIDs.begin(); itr != colliIDs.end(); ++itr) {
		if (itr->second->size() <= opt::maxGroupSize)
			m_colliIDs[itr->first] = itr->second;
		else
			m_colliIDs[itr->first] = NULL;
	}
	idFH.close();

	m_countTable = vector<size_t>(m_fullIDs.size(), 0.0);

	string countsFilename = (filterFile).substr(0, (filterFile).length() - 3)
			+ "_indexCount.tsv";

	//load in freqTable of seeds
	ifstream freqFH(countsFilename.c_str(), ifstream::in);
	getline(freqFH, line);
	while (freqFH.good()) {
		ID id;
		stringstream converter(line);
		converter >> id;
		size_t count;
		converter >> count;
		m_countTable[id] = count;
		getline(freqFH, line);
	}
	freqFH.close();

	size_t sum = 0;

	for (vector<size_t>::const_iterator itr = m_countTable.begin();
			itr != m_countTable.end(); ++itr) {
		sum += *itr;
	}
	m_freqTable = vector<double>(m_fullIDs.size(), 0.0);
	m_perFrameProb = vector<double>(m_fullIDs.size(), 0.0);
	for (size_t i = 0;
			i < m_countTable.size(); ++i) {
		m_freqTable[i] = double(m_countTable[i])/double(sum);
		m_perFrameProb[i] = calcProbSingleFrame(m_freqTable[i]);
//		cerr << i << "\t" << m_freqTable[i] << "\t" << m_perFrameProb[i] << "\t"
//				<< getMinCount(91, m_perFrameProb[i]) << endl;
	}
	if(opt::allowMisses > 0 && m_filter.getSeedValues().size()  == 0)
	{
		cerr << "Allowed miss (-a) should not be used with k-mers" << endl;
		exit(1);
	}
	m_baseProb = m_filter.getFPR(opt::allowMisses);
}

void BloomMapClassifier::filter(const vector<string> &inputFiles) {
	string outputName = opt::outputPrefix + "_reads.tsv";

	//TODO output fasta?
	if (opt::outputType == "fq") {
		outputName = opt::outputPrefix + "_reads.fq";
	}

	Dynamicofstream readsOutput(outputName);

	//print out header info and initialize variables
	ResultsManager<ID> resSummary(m_fullIDs, false);

	cerr << "Filtering Start" << endl;

	//debugging
	size_t multiCount = 0;
	size_t incorrectCount = 0;

	for (vector<string>::const_iterator it = inputFiles.begin();
			it != inputFiles.end(); ++it) {
		gzFile fp;
		fp = gzopen(it->c_str(), "r");
		if (fp == NULL) {
			cerr << "Cannot open file" << it->c_str() << endl;
			exit(1);
		}
		kseq_t *seq = kseq_init(fp);
		int l;
#pragma omp parallel private(l)
		for (string sequence;;) {
			string name;
			string qual;
			string comment;
#pragma omp critical(sequence)
			{
				l = kseq_read(seq);
				sequence = string(seq->seq.s, seq->seq.l);
				name = string(seq->name.s, seq->name.l);
				qual = string(seq->qual.s, seq->qual.l);
				comment = string(seq->comment.s, seq->comment.l);
			}
			if (l >= 0) {
#pragma omp critical(totalReads)
				if (++m_numRead % opt::fileInterval == 0) {
					cerr << "Currently Reading Read Number: " << m_numRead
							<< endl;
				}
				vector<vector<ID> > hitsPattern;
				vector<unsigned> saturation;
				unsigned evaluatedSeeds = 0;
				vector<unsigned> sig = getMatchSignature(sequence,
						evaluatedSeeds, hitsPattern, saturation);
				double pVal = 1.0;
				unsigned maxCount = 0;
				vector<ID> signifResults;
				vector<unsigned> signifValues;
				vector<unsigned> fullSignifCounts;
				ID idIndex = evalRead(hitsPattern, evaluatedSeeds, pVal,
						maxCount, signifResults, signifValues, fullSignifCounts);
				resSummary.updateSummaryData(idIndex);
				const string &fullID =
						idIndex == opt::EMPTY ? NO_MATCH :
						idIndex == resSummary.getMultiMatchIndex() ?
								MULTI_MATCH : m_fullIDs.at(idIndex);
				//debugging
//#pragma omp critical(cout)
//				if (name == "noMatch") {
//					cout << pVal << "\t" << "F" << "\t" << name << "\t"
//							<< fullID << "\t" << probMulti << endl;
//				} else {
//					cout << pVal << "\t" << "T" << "\t" << name << "\t"
//							<< fullID << "\t" << probMulti << endl;
//				}
				bool match = false;
				unsigned tempCount = 0;
				for (unsigned i = 0; i < signifResults.size(); ++i) {
					//TODO resolve precision error better
					if (fullSignifCounts[i] + 3 >= maxCount) {
						if (name == m_fullIDs.at(signifResults[i])) {
							match = true;
						}
						++tempCount;
					}
				}
				if(tempCount > 1){
#pragma omp atomic
					++multiCount;
				}

#pragma omp critical(cout)
				if (!match) {
#pragma omp atomic
					++incorrectCount;
//				if (name != fullID) {
//					unsigned correctCount = 0;
//					for (unsigned i = 0; i < signifResults.size(); ++i) {
//						if(  signifResults[i] == m_idToIndex[name]){
//							correctCount = signifValues[i];
//							break;
//						}
//					}
//					cout << maxCount - correctCount << endl;
					cout << m_numRead << "\tCorrectID:" << m_idToIndex[name]
							<< "\tCorrectName:" << name << "\tPredictedID:"
							<< m_idToIndex[fullID] << "\tPredictedName:"
							<< fullID << "\tCorrectID:"
							<< base64_chars[m_idToIndex[name] % 64] << "\tpVal:"
							<< log10(pVal) * (-10.0) << "\t" << endl;

					for (unsigned i = 0; i < signifResults.size(); ++i) {
						cout << signifResults[i] << ","
								<< m_fullIDs[signifResults[i]] << ","
								<< base64_chars[signifResults[i] % 64] << ","
								<< signifValues[i] << "," << fullSignifCounts[i]
								<< "," << m_freqTable[signifResults[i]] << "\t";
					}
					cout << endl;
					printVerbose(name, comment, sequence, hitsPattern, sig,
							saturation, idIndex);
				}

				if (idIndex != opt::EMPTY) {
					if (opt::outputType == "fq") {
#pragma omp critical(outputFiles)
						{
							readsOutput << "@" << name << " " << fullID << "\n"
									<< sequence << "\n+\n" << qual << "\n";
						}
					} else {
#pragma omp critical(outputFiles)
						{
							readsOutput << fullID << "\t" << name << "\t";
							readsOutput << "\n";
						}
					}
				}
			} else {
				break;
			}
		}
		kseq_destroy(seq);
		gzclose(fp);
	}

	cerr << "Multiple Map Count:" << multiCount << endl;
	cerr << incorrectCount << endl;

	cerr << "Total Reads:" << m_numRead << endl;
	cerr << "Writing file: " << opt::outputPrefix.c_str() << "_summary.tsv"
			<< endl;

	Dynamicofstream summaryOutput(opt::outputPrefix + "_summary.tsv");
	summaryOutput << resSummary.getResultsSummary(m_numRead);
	summaryOutput.close();
	cout.flush();
}

void BloomMapClassifier::filterPair(const string &file1, const string &file2) {
	string outputName = opt::outputPrefix + "_reads.tsv";

	//TODO output fasta?
	if (opt::outputType == "fq") {
		outputName = opt::outputPrefix + "_reads.fq";
	}

	Dynamicofstream readsOutput(outputName);

	//print out header info and initialize variables
	ResultsManager<ID> resSummary(m_fullIDs, false);

	if(!m_filter.getSeedValues().empty()){
		cerr << "Spaced Seed Mode" << endl;
	}

	cerr << "Filtering Start" << endl;

	gzFile fp1;
	gzFile fp2;
	fp1 = gzopen(file1.c_str(), "r");
	if (fp1 == NULL) {
		cerr << "Cannot open file " << file1.c_str() << endl;
		exit(1);
	}
	fp2 = gzopen(file2.c_str(), "r");
	if (fp2 == NULL) {
		cerr << "Cannot open file " << file2.c_str() << endl;
		exit(1);
	}
	kseq_t *seq1 = kseq_init(fp1);
	kseq_t *seq2 = kseq_init(fp2);
	int l1;
	int l2;
#pragma omp parallel private(l1, l2)
	for (string sequence1;;) {
		string sequence2;
		//TODO sanity check if names to make sure both names are the same
		string name1;
		string name2;
		string qual1;
		string qual2;
#pragma omp critical(sequence)
		{
			l1 = kseq_read(seq1);
			sequence1 = string(seq1->seq.s, seq1->seq.l);
			name1 = string(seq1->name.s, seq1->name.l);
			qual1 = string(seq1->qual.s, seq1->qual.l);

			l2 = kseq_read(seq2);
			sequence2 = string(seq2->seq.s, seq2->seq.l);
			name2 = string(seq2->name.s, seq2->name.l);
			qual2 = string(seq2->qual.s, seq2->qual.l);
		}
		if (l1 >= 0 && l2 >= 0) {
#pragma omp critical(totalReads)
			if (++m_numRead % opt::fileInterval == 0) {
				cerr << "Currently Reading Read Number: " << m_numRead << endl;
			}

			google::dense_hash_map<ID, unsigned> hitCounts1;
			google::dense_hash_map<ID, unsigned> hitCounts2;
			hitCounts1.set_empty_key(opt::EMPTY);
			hitCounts2.set_empty_key(opt::EMPTY);

//			unsigned evaluatedSeeds = 0;
//			unsigned rmCount = 0;
//			double probFP = 0.0;
//			vector<vector<ID> > hitsPattern;
//			vector<unsigned> sig = getMatchSignature(sequence1, evaluatedSeeds,
//					hitsPattern);
//			vector<unsigned> rmMatch = calcProb(sig, rmCount, evaluatedSeeds,
//					probFP);
//			if (rmMatch.size() > 0) {
//				printVerbose(name1, sequence1, sig, evaluatedSeeds, rmCount,
//						rmMatch, probFP, hitsPattern);
//			}

			if(m_numRead > opt::fileInterval){
				exit(1);
			}

//			unsigned score1 = evaluateRead(sequence1, hitCounts1);
//			unsigned score2 = evaluateRead(sequence2, hitCounts2);
			unsigned bestHit = 0;

			vector<ID> hits;
//			if (opt::inclusive) {
//				if (score1 > threshold1 || score2 > threshold2)
//					bestHit = convertToHitsOnlyOne(hitCounts1, hitCounts2,
//							hits);
//			} else {
//				if (score1 > threshold1 && score2 > threshold2)
//					bestHit = convertToHitsBoth(hitCounts1, hitCounts2, hits);
//			}

			ID idIndex = opt::EMPTY;

			idIndex = resSummary.updateSummaryData(hits);
			if (idIndex != opt::EMPTY) {
				const string &fullID =
						idIndex == opt::EMPTY ? NO_MATCH :
						idIndex == resSummary.getMultiMatchIndex() ?
								MULTI_MATCH : m_fullIDs.at(idIndex);
				if (opt::outputType == "fq") {
#pragma omp critical(outputFiles)
					{
						readsOutput << "@" << name1 << " " << fullID << "\n"
								<< sequence1 << "\n+\n" << qual1 << "\n" << "@"
								<< name2 << " " << fullID << "\n" << sequence2
								<< "\n+\n" << qual2 << "\n";
					}
				} else {
#pragma omp critical(outputFiles)
					{
						readsOutput << fullID << "\t" << name1 << "\t"
								<< bestHit << "\n";
					}
				}
			}
		} else {
			if (l1 != -1 || l2 != -1) {
				cerr << "Terminated without getting to eof at read "
						<< m_numRead << endl;
				exit(1);
			}
			break;
		}
	}
	kseq_destroy(seq1);
	kseq_destroy(seq2);
	gzclose(fp1);
	gzclose(fp2);

	cerr << "Total Reads:" << m_numRead << endl;
	cerr << "Writing file: " << opt::outputPrefix.c_str() << "_summary.tsv"
			<< endl;

	Dynamicofstream summaryOutput(opt::outputPrefix + "_summary.tsv");
	summaryOutput << resSummary.getResultsSummary(m_numRead);
	summaryOutput.close();
	cout.flush();
}

BloomMapClassifier::~BloomMapClassifier() {
}

