/*
 * BloomMapClassifier.cpp
 *
 *  Created on: Mar 23, 2016
 *      Author: cjustin
 */

#include <BloomMapClassifier.h>
#include "DataLayer/FastaReader.h"
#include "Common/Options.h"
#include "Options.h"
#include <iostream>
//#include "Common/Dynamicofstream.h"
#include "ResultsManager.h"
#include <zlib.h>
#include "DataLayer/kseq.h"
KSEQ_INIT(gzFile, gzread)

BloomMapClassifier::BloomMapClassifier(const string &filterFile) :
		m_filter(BloomMapSSBitVec<ID>(filterFile)) {
	cerr << "FPR given allowed misses: "
			<< m_filter.getFPR(opt::allowMisses)
			<< endl;
	//load in ID file
	string idFile = (filterFile).substr(0, (filterFile).length() - 3)
			+ "_ids.txt";
	google::dense_hash_map<ID, string> ids;
	google::dense_hash_map<ID, boost::shared_ptr< vector<ID> > > colliIDs;
	colliIDs.set_empty_key(opt::EMPTY);
	ids.set_empty_key(opt::EMPTY);
	if (!fexists(idFile)) {
		cerr << "Error: " << idFile.c_str()
				<< " File cannot be opened. A corresponding id file is needed."
				<< endl;
		exit(1);
	}
	else{
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
		for (;converter >> fullID; ++idCount) {
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
	m_fullIDs[0] = "missing"; //set first element to missing
	for (google::dense_hash_map<ID, string>::const_iterator itr = ids.begin();
			itr != ids.end(); ++itr) {
		m_fullIDs[itr->first] = itr->second;
	}
	m_colliIDs = vector< boost::shared_ptr<vector<ID> > >(colliIDs.size() + 1);
	for (google::dense_hash_map<ID, boost::shared_ptr< vector<ID> > >::const_iterator itr =
			colliIDs.begin(); itr != colliIDs.end(); ++itr) {
		m_colliIDs[itr->first] = itr->second;
	}
	idFH.close();
}

void BloomMapClassifier::filter(const vector<string> &inputFiles) {
	size_t totalReads = 0;

	Dynamicofstream readsOutput(opt::outputPrefix + "_reads.tsv");

	//print out header info and initialize variables
	ResultsManager resSummary(m_fullIDs, false);

	cerr << "Filtering Start" << endl;

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
#pragma omp critical(sequence)
			{
				l = kseq_read(seq);
				sequence = string(seq->seq.s, seq->seq.l);
				name = string(seq->name.s, seq->name.l);
				qual = string(seq->qual.s, seq->qual.l);
			}
			if (l >= 0) {
#pragma omp atomic
				++totalReads;
#pragma omp critical(totalReads)
				if (totalReads % 1000000 == 0) {
					cerr << "Currently Reading Read Number: " << totalReads
							<< endl;
				}

				google::dense_hash_map<ID, unsigned> hitCounts;
				hitCounts.set_empty_key(opt::EMPTY);
				google::dense_hash_map<ID, unsigned> hitCountsSolid;
				hitCountsSolid.set_empty_key(opt::EMPTY);
				unsigned solidCount1 = 0;
				unsigned score = evaluateRead(name, hitCounts, hitCountsSolid, solidCount1);
				assert(score);
				unsigned threshold = opt::score
						* (l - m_filter.getKmerSize() + 1);
				vector<ID> hits;
				bool aboveThreshold = score > threshold;
				if (aboveThreshold) {
					convertToHits(hitCounts, hits);
				}
				ID idIndex = resSummary.updateSummaryData(hits, aboveThreshold);
				if (idIndex != opt::EMPTY) {
					const string &fullID =
							idIndex == opt::COLLI ?
									UNKNOWN : m_fullIDs.at(idIndex);
					if (opt::outputType == "fq") {
#pragma omp critical(outputFiles)
						{
							readsOutput << "@" << name << " " << fullID << "\n"
									<< sequence << "\n+\n" << qual << "\n";
						}
					} else {
#pragma omp critical(outputFiles)
						{
							readsOutput << fullID << "\t" << name << "\n";
						}
					}
				}
			} else {
				kseq_destroy(seq);
				break;
			}
		}
		gzclose(fp);
	}

	cerr << "Total Reads:" << totalReads << endl;
	cerr << "Writing file: " << opt::outputPrefix.c_str() << "_summary.tsv" << endl;

	Dynamicofstream summaryOutput(opt::outputPrefix + "_summary.tsv");
	summaryOutput << resSummary.getResultsSummary(totalReads);
	summaryOutput.close();
	cout.flush();
}

void BloomMapClassifier::filterPair(const string &file1, const string &file2) {
	size_t totalReads = 0;
	string outputName = opt::outputPrefix + "_reads.tsv";

	//TODO output fasta?
	if(opt::outputType == "fq"){
		outputName = opt::outputPrefix + "_reads.fq";
	}

	Dynamicofstream readsOutput(outputName);

	//print out header info and initialize variables
	ResultsManager resSummary(m_fullIDs, false);

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
#pragma omp atomic
			++totalReads;
#pragma omp critical(totalReads)
			if (totalReads % 1000000 == 0) {
				cerr << "Currently Reading Read Number: " << totalReads << endl;
			}

			google::dense_hash_map<ID, unsigned> hitCounts1;
			google::dense_hash_map<ID, unsigned> hitCounts2;
//			google::dense_hash_map<ID, unsigned> hitCountsSolid1;
//			google::dense_hash_map<ID, unsigned> hitCountsSolid2;
			hitCounts1.set_empty_key(opt::EMPTY);
			hitCounts2.set_empty_key(opt::EMPTY);
//			hitCountsSolid1.set_empty_key(opt::EMPTY);
//			hitCountsSolid2.set_empty_key(opt::EMPTY);

			unsigned solidCount1 = 0;
			unsigned solidCount2 = 0;

//			unsigned score1 = evaluateRead(sequence1, hitCounts1,
//					hitCountsSolid1, solidCount1);
//			unsigned score2 = evaluateRead(sequence2, hitCounts2,
//					hitCountsSolid2, solidCount2);
			unsigned score1 = evaluateRead(sequence1, hitCounts1);
			unsigned score2 = evaluateRead(sequence2, hitCounts2);
			unsigned threshold1 =
					opt::minHitOnly ?
							opt::minHit - 1 :
							opt::score * (l1 - m_filter.getKmerSize() + 1);
			unsigned threshold2 = opt::score
					* (l2 - m_filter.getKmerSize() + 1);

			vector<ID> hits;
			vector<ID> solidHits;
			bool aboveThreshold = false;
			unsigned bestScoreDiff = 0;
			unsigned bestScoreDiffSolid = 0;
			if (opt::inclusive) {
				aboveThreshold = score1 >= threshold1 || score2 >= threshold2;
				if (aboveThreshold) {
					convertToHitsOnlyOne(hitCounts1, hitCounts2, hits,
							bestScoreDiff);
//					convertToHitsOnlyOne(hitCountsSolid1, hitCountsSolid2,
//							solidHits, bestScoreDiffSolid);
				}
			} else {
				aboveThreshold = score1 >= threshold1 && score2 >= threshold2;
				if (aboveThreshold) {
					convertToHitsBoth(hitCounts1, hitCounts2, hits);
//					convertToHitsBoth(hitCountsSolid1, hitCountsSolid2,
//							solidHits);
				}
			}

			ID idIndex = 0;
			//compare solid hits with hits
			//if solidHits are the same or have no hits return results
//			if (hits.size() == solidHits.size() || solidCount1 + solidCount2 == 0 || solidHits.size() == 0) {

			//TODO: make into option!
//			unsigned delta2 = 1;


			if (bestScoreDiff > opt::delta) {
				idIndex = resSummary.updateSummaryData(hits, aboveThreshold);
			} else {
				hits.clear();
				idIndex = resSummary.updateSummaryData(hits, aboveThreshold);
//				if (bestScoreDiffSolid <= delta2) {
//					solidHits.clear();
//				}
//				idIndex = resSummary.updateSummaryData(solidHits,
//						aboveThreshold);
			}
			//if solidHits differ, assign ID to solidHit IDs (even if more ambiguous mapping)

//				if (idIndex != opt::EMPTY) {
			const string &fullID =
					idIndex == opt::COLLI ? UNKNOWN : m_fullIDs.at(idIndex);
			if (opt::outputType == "fq") {
#pragma omp critical(outputFiles)
				{
					readsOutput << "@" << name1 << " " << fullID << "\n"
							<< sequence1 << "\n+\n" << qual1 << "\n";
					readsOutput << "@" << name2 << " " << fullID << "\n"
							<< sequence2 << "\n+\n" << qual2 << "\n";
				}
			} else {
#pragma omp critical(outputFiles)
				{
					readsOutput << fullID << "\t" << name1 << "\t" << score1
							<< "\t" << score2 << "\t" << solidCount1 << "\t"
							<< solidCount2 << "\t" << bestScoreDiff << "\t" << bestScoreDiffSolid;
//					google::dense_hash_set<ID> tempSet;
//					tempSet.set_empty_key(opt::EMPTY);
//					for (google::dense_hash_map<ID, unsigned>::const_iterator i =
//							hitCounts1.begin(); i != hitCounts1.end(); ++i) {
//						google::dense_hash_map<ID, unsigned>::const_iterator itr =
//								hitCounts2.find(i->first);
//						if (itr != hitCounts2.end()) {
//							tempSet.insert(i->first);
//							readsOutput << "\t" << i->first << ":" << i->second
//									<< ";" << itr->second;
//						}
//						else {
//							readsOutput << "\t" << i->first << ":" << i->second
//									<< ";0";
//						}
//					}
//					for (google::dense_hash_map<ID, unsigned>::const_iterator i =
//							hitCounts2.begin(); i != hitCounts2.end(); ++i) {
//						if (tempSet.find(i->first) == tempSet.end()) {
//							readsOutput << "\t" << i->first << ":0;"
//									<< i->second;
//						}
//					}
//					readsOutput << "\t|";
//					google::dense_hash_set<ID> tempSet2;
//					tempSet2.set_empty_key(opt::EMPTY);
//					for (google::dense_hash_map<ID, unsigned>::const_iterator i =
//							hitCountsSolid1.begin(); i != hitCountsSolid1.end();
//							++i) {
//						google::dense_hash_map<ID, unsigned>::const_iterator itr =
//								hitCountsSolid2.find(i->first);
//						if (itr != hitCountsSolid2.end()) {
//							tempSet2.insert(i->first);
//							readsOutput << "\t" << i->first << ":" << i->second
//									<< ";" << itr->second;
//						} else {
//							readsOutput << "\t" << i->first << ":" << i->second
//									<< ";0";
//						}
//					}
//					for (google::dense_hash_map<ID, unsigned>::const_iterator i =
//							hitCountsSolid2.begin(); i != hitCountsSolid2.end();
//							++i) {
//						if (tempSet2.find(i->first) == tempSet2.end()) {
//							readsOutput << "\t" << i->first << ":0;"
//									<< i->second;
//						}
//					}
					readsOutput << "\n";
				}
//					}
			}
		} else {
			if (l1 != -1 || l2 != -1) {
				cerr << "Terminated without getting to eof at read "
						<< totalReads << endl;
				exit(1);
			}
			break;
		}
	}
	kseq_destroy(seq1);
	kseq_destroy(seq2);
	gzclose(fp1);
	gzclose(fp2);

	cerr << "Total Reads:" << totalReads << endl;
	cerr << "Writing file: " << opt::outputPrefix.c_str() << "_summary.tsv" << endl;

	Dynamicofstream summaryOutput(opt::outputPrefix + "_summary.tsv");
	summaryOutput << resSummary.getResultsSummary(totalReads);
	summaryOutput.close();
	cout.flush();
}

BloomMapClassifier::~BloomMapClassifier() {
}

