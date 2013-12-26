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
#include "Common/Dynamicofstream.h"
#include "ResultsManager.h"

BioBloomClassifier::BioBloomClassifier(const vector<string> &filterFilePaths,
		size_t minHit, double percentMinHit, size_t maxHitValue,
		const string &prefix, const string &outputPostFix, uint8_t tileModifier) :
		minHit(minHit), percentMinHit(percentMinHit), filterNum(
				filterFilePaths.size()), maxHitValue(maxHitValue), noMatch(
				"noMatch"), multiMatch("multiMatch"), prefix(prefix), postfix(
				outputPostFix), tileModifier(tileModifier) {
	loadFilters(filterFilePaths);
}

/*
 * Generic filtering function (signel end, no fa or fq file outputs)
 */
void BioBloomClassifier::filter(const vector<string> &inputFiles) {
	Dynamicofstream readStatusOutput(prefix + "_status.tsv" + postfix);

	//results summary object
	ResultsManager resSummary(hashSigs, filters, infoFiles, minHit,
			percentMinHit, maxHitValue, tileModifier);

	size_t totalReads = 0;

	//print out header info and initialize variables
	readStatusOutput << getReadSummaryHeader(hashSigs);

	cerr << "Filtering Start" << endl;

	if (tileModifier == 0) {
		for (vector<string>::const_iterator it = inputFiles.begin();
				it != inputFiles.end(); ++it) {
			FastaReader sequence(it->c_str(), FastaReader::NO_FOLD_CASE);
			FastqRecord rec;

			//stored out of loop so reallocation does not have to be done
			unordered_map<string, size_t> hits(filterNum);
			while (sequence >> rec) {

				//track read progress
				++totalReads;
				if (totalReads % 1000000 == 0) {
					cerr << "Currently Reading Read Number: " << totalReads
							<< endl;
				}

				//initialize hits to zero
				initHits(hits);

				//for each hashSigniture/kmer combo multi, cut up read into kmer sized used
				for (vector<string>::const_iterator j = hashSigs.begin();
						j != hashSigs.end(); ++j) {
					evaluateRead(rec, *j, hits);
				}

				//print hit results to read status
				readStatusOutput
						<< getReadStatStr(rec.id, rec.seq.length(), hits);

				//Evaluate hit data and record for summary
				resSummary.updateSummaryData(rec.seq.length(), hits);
			}
		}
	} else {
		for (vector<string>::const_iterator it = inputFiles.begin();
				it != inputFiles.end(); ++it) {
			FastaReader sequence(it->c_str(), FastaReader::NO_FOLD_CASE);
			FastqRecord rec;

			//stored out of loop so reallocation does not have to be done
			unordered_map<string, size_t> hits(filterNum);
			while (sequence >> rec) {

				//track read progress
				++totalReads;
				if (totalReads % 1000000 == 0) {
					cerr << "Currently Reading Read Number: " << totalReads
							<< endl;
				}

				//initialize hits to zero
				initHits(hits);

				//for each hashSigniture/kmer combo multi, cut up read into kmer sized used
				for (vector<string>::const_iterator j = hashSigs.begin();
						j != hashSigs.end(); ++j) {
					evaluateRead(rec, *j, hits, tileModifier);
				}

				//print hit results to read status
				readStatusOutput
						<< getReadStatStr(rec.id, rec.seq.length(), hits);

				//Evaluate hit data and record for summary
				resSummary.updateSummaryData(rec.seq.length(), hits);
			}
		}
	}

	readStatusOutput.close();

	cerr << "Total Reads:" << totalReads << endl;

	Dynamicofstream summaryOutput(prefix + "_summary.tsv");
	summaryOutput << resSummary.getResultsSummary(totalReads);
	summaryOutput.close();

	if (maxHitValue > 0) {
		Dynamicofstream countSummaryOutput(prefix + "_rawCounts.tsv");
		countSummaryOutput << resSummary.getCountSummary(totalReads);
		countSummaryOutput.close();
	}
}

/*
 * Filters reads
 * Assumes only one hash signature exists (load only filters with same
 * hash functions)
 * Prints reads into seperate files
 *
 * outputType must be fa or fq
 */
void BioBloomClassifier::filterPrint(const vector<string> &inputFiles,
		const string &outputType) {

	Dynamicofstream readStatusOutput(prefix + "_status.tsv" + postfix);

	//results summary object
	ResultsManager resSummary(hashSigs, filters, infoFiles, minHit,
			percentMinHit, maxHitValue, tileModifier);

	size_t totalReads = 0;

	unordered_map<string, shared_ptr<Dynamicofstream> > outputFiles;
	shared_ptr<Dynamicofstream> no_match(
			new Dynamicofstream(
					prefix + "_" + noMatch + "." + outputType + postfix));
	shared_ptr<Dynamicofstream> multi_match(
			new Dynamicofstream(
					prefix + "_" + multiMatch + "." + outputType + postfix));
	outputFiles[noMatch] = no_match;
	outputFiles[multiMatch] = multi_match;

	//initialize variables
	for (vector<string>::const_iterator j = hashSigs.begin();
			j != hashSigs.end(); ++j) {
		const vector<string> idsInFilter = (*filters[*j]).getFilterIds();
		for (vector<string>::const_iterator i = idsInFilter.begin();
				i != idsInFilter.end(); ++i) {
			shared_ptr<Dynamicofstream> temp(
					new Dynamicofstream(
							prefix + "_" + *i + "." + outputType + postfix));
			outputFiles[*i] = temp;
		}
	}

	//print out header info and initialize variables
	readStatusOutput << getReadSummaryHeader(hashSigs);

	cerr << "Filtering Start" << endl;

	if (tileModifier == 0) {

		for (vector<string>::const_iterator it = inputFiles.begin();
				it != inputFiles.end(); ++it) {
			FastaReader sequence(it->c_str(), FastaReader::NO_FOLD_CASE);
			FastqRecord rec;
			//hits results stored in hashmap of filternames and hits
			unordered_map<string, size_t> hits(filterNum);
			while (sequence >> rec) {

				++totalReads;
				if (totalReads % 1000000 == 0) {
					cerr << "Currently Reading Read Number: " << totalReads
							<< endl;
				}

				//initialize hits to zero
				initHits(hits);

				//for each hashSigniture/kmer combo multi, cut up read into kmer sized used
				for (vector<string>::const_iterator j = hashSigs.begin();
						j != hashSigs.end(); ++j) {
					evaluateRead(rec, *j, hits);
				}

				//print hit results to read status
				readStatusOutput
						<< getReadStatStr(rec.id, rec.seq.length(), hits);

				//Evaluate hit data and record for summary
				const string &outputFileName = resSummary.updateSummaryData(
						rec.seq.length(), hits);

				if (outputType == "fa") {
					(*outputFiles[outputFileName]) << ">" << rec.id << "\n"
							<< rec.seq << "\n";
				} else {
					(*outputFiles[outputFileName]) << "@" << rec.id << "\n"
							<< rec.seq << "\n+\n" << rec.qual << "\n";
				}
			}
		}
	} else {

		for (vector<string>::const_iterator it = inputFiles.begin();
				it != inputFiles.end(); ++it) {
			FastaReader sequence(it->c_str(), FastaReader::NO_FOLD_CASE);
			FastqRecord rec;
			//hits results stored in hashmap of filternames and hits
			unordered_map<string, size_t> hits(filterNum);
			while (sequence >> rec) {
				++totalReads;
				if (totalReads % 1000000 == 0) {
					cerr << "Currently Reading Read Number: " << totalReads
							<< endl;
				}

				//initialize hits to zero
				initHits(hits);

				//for each hashSigniture/kmer combo multi, cut up read into kmer sized used
				for (vector<string>::const_iterator j = hashSigs.begin();
						j != hashSigs.end(); ++j) {
					evaluateRead(rec, *j, hits, tileModifier);
				}

				//print hit results to read status
				readStatusOutput
						<< getReadStatStr(rec.id, rec.seq.length(), hits);

				//Evaluate hit data and record for summary
				const string &outputFileName = resSummary.updateSummaryData(
						rec.seq.length(), hits);

				if (outputType == "fa") {
					(*outputFiles[outputFileName]) << ">" << rec.id << "\n"
							<< rec.seq << "\n";
				} else {
					(*outputFiles[outputFileName]) << "@" << rec.id << "\n"
							<< rec.seq << "\n+\n" << rec.qual << "\n";
				}
			}
		}
	}

	readStatusOutput.close();

	//close sorting files
	for (unordered_map<string, shared_ptr<Dynamicofstream> >::iterator j =
			outputFiles.begin(); j != outputFiles.end(); ++j) {
		j->second->close();
	}
	cerr << "Total Reads:" << totalReads << endl;

	Dynamicofstream summaryOutput(prefix + "_summary.tsv");
	summaryOutput << resSummary.getResultsSummary(totalReads);
	summaryOutput.close();

	if (maxHitValue > 0) {
		Dynamicofstream countSummaryOutput(prefix + "_rawCounts.tsv");
		countSummaryOutput << resSummary.getCountSummary(totalReads);
		countSummaryOutput.close();
	}
}

/*
 * Filters reads -> uses paired end information
 * Assumes only one hash signature exists (load only filters with same
 * hash functions)
 */
void BioBloomClassifier::filterPair(const string &file1, const string &file2) {

	Dynamicofstream readStatusOutput(prefix + "_status.tsv" + postfix);

	//results summary object
	ResultsManager resSummary(hashSigs, filters, infoFiles, minHit,
			percentMinHit, maxHitValue, tileModifier);

	size_t totalReads = 0;

	//print out header info and initialize variables for summary
	readStatusOutput << getReadSummaryHeader(hashSigs);

	cerr << "Filtering Start" << "\n";

	FastaReader sequence1(file1.c_str(), FastaReader::NO_FOLD_CASE);
	FastaReader sequence2(file2.c_str(), FastaReader::NO_FOLD_CASE);
	FastqRecord rec1;
	FastqRecord rec2;
	//hits results stored in hashmap of filter names and hits
	unordered_map<string, size_t> hits1(filterNum);
	unordered_map<string, size_t> hits2(filterNum);

	if (tileModifier == 0) {
		while (sequence1 >> rec1 && sequence2 >> rec2) {
			++totalReads;
			if (totalReads % 1000000 == 0) {
				cerr << "Currently Reading Read Number: " << totalReads << endl;
			}

			//initialize hits to zero
			initHits(hits1);
			initHits(hits2);

			//for each hashSigniture/kmer combo multi, cut up read into kmer sized used
			for (vector<string>::const_iterator j = hashSigs.begin();
					j != hashSigs.end(); ++j) {
				string tempStr1 = rec1.id.substr(0, rec1.id.find_last_of("/"));
				string tempStr2 = rec2.id.substr(0, rec2.id.find_last_of("/"));
				if (tempStr1 == tempStr2) {
					evaluateRead(rec1, *j, hits1);
					evaluateRead(rec2, *j, hits2);
				} else {
					cerr << "Read IDs do not match" << "\n" << tempStr1 << "\n"
							<< tempStr2 << endl;
					exit(1);
				}
			}

			string readID = rec1.id.substr(0, rec1.id.length() - 2);

			//print hit results to read status
			readStatusOutput
					<< getReadStatStrPair(readID, rec1.seq.length(),
							rec2.seq.length(), hits1, hits2);

			//Evaluate hit data and record for summary
			resSummary.updateSummaryData(rec1.seq.length(), rec2.seq.length(),
					hits1, hits2);
		}
	} else {

		while (sequence1 >> rec1 && sequence2 >> rec2) {
			++totalReads;
			if (totalReads % 1000000 == 0) {
				cerr << "Currently Reading Read Number: " << totalReads << endl;
			}

			//initialize hits to zero
			initHits(hits1);
			initHits(hits2);

			//for each hashSigniture/kmer combo multi, cut up read into kmer sized used
			for (vector<string>::const_iterator j = hashSigs.begin();
					j != hashSigs.end(); ++j) {
				string tempStr1 = rec1.id.substr(0, rec1.id.find_last_of("/"));
				string tempStr2 = rec2.id.substr(0, rec2.id.find_last_of("/"));
				if (tempStr1 == tempStr2) {
					evaluateRead(rec1, *j, hits1, tileModifier);
					evaluateRead(rec2, *j, hits2, tileModifier);
				} else {
					cerr << "Read IDs do not match" << "\n" << tempStr1 << "\n"
							<< tempStr2 << endl;
					exit(1);
				}
			}

			string readID = rec1.id.substr(0, rec1.id.length() - 2);

			//print hit results to read status
			readStatusOutput
					<< getReadStatStrPair(readID, rec1.seq.length(),
							rec2.seq.length(), hits1, hits2);

			//Evaluate hit data and record for summary
			resSummary.updateSummaryData(rec1.seq.length(), rec2.seq.length(),
					hits1, hits2);
		}
	}

	if (sequence2 >> rec2 && sequence1.eof() && sequence2.eof()) {
		cerr
				<< "error: eof bit not flipped. Input files may be different lengths"
				<< endl;
	}

	readStatusOutput.close();

	cerr << "Total Reads:" << totalReads << endl;

	Dynamicofstream summaryOutput(prefix + "_summary.tsv");
	summaryOutput << resSummary.getResultsSummary(totalReads);
	summaryOutput.close();

	if (maxHitValue > 0) {
		Dynamicofstream countSummaryOutput(prefix + "_rawCounts.tsv");
		countSummaryOutput << resSummary.getCountSummary(totalReads);
		countSummaryOutput.close();
	}
}

/*
 * Filters reads -> uses paired end information
 * Assumes only one hash signature exists (load only filters with same
 * hash functions)
 * prints reads
 */
void BioBloomClassifier::filterPairPrint(const string &file1,
		const string &file2, const string &outputType) {
	Dynamicofstream readStatusOutput(prefix + "_status.tsv" + postfix);

	//results summary object
	ResultsManager resSummary(hashSigs, filters, infoFiles, minHit,
			percentMinHit, maxHitValue, tileModifier);

	size_t totalReads = 0;

	unordered_map<string, shared_ptr<Dynamicofstream> > outputFiles;
	shared_ptr<Dynamicofstream> noMatch1(
			new Dynamicofstream(
					prefix + "_" + noMatch + "_1." + outputType + postfix));
	shared_ptr<Dynamicofstream> noMatch2(
			new Dynamicofstream(
					prefix + "_" + noMatch + "_2." + outputType + postfix));
	shared_ptr<Dynamicofstream> multiMatch1(
			new Dynamicofstream(
					prefix + "_" + multiMatch + "_1." + outputType + postfix));
	shared_ptr<Dynamicofstream> multiMatch2(
			new Dynamicofstream(
					prefix + "_" + multiMatch + "_2." + outputType + postfix));
	outputFiles[noMatch + "1"] = noMatch1;
	outputFiles[noMatch + "2"] = noMatch2;
	outputFiles[multiMatch + "1"] = multiMatch1;
	outputFiles[multiMatch + "2"] = multiMatch2;

	//initialize variables and print filter ids
	for (vector<string>::const_iterator j = hashSigs.begin();
			j != hashSigs.end(); ++j) {
		const vector<string> idsInFilter = (*filters[*j]).getFilterIds();
		for (vector<string>::const_iterator i = idsInFilter.begin();
				i != idsInFilter.end(); ++i) {
			shared_ptr<Dynamicofstream> temp1(
					new Dynamicofstream(
							prefix + "_" + *i + "_1." + outputType + postfix));
			shared_ptr<Dynamicofstream> temp2(
					new Dynamicofstream(
							prefix + "_" + *i + "_2." + outputType + postfix));
			outputFiles[*i + "1"] = temp1;
			outputFiles[*i + "2"] = temp2;
		}
	}

	//print out header info and initialize variables for summary
	readStatusOutput << getReadSummaryHeader(hashSigs);

	cerr << "Filtering Start" << "\n";

	FastaReader sequence1(file1.c_str(), FastaReader::NO_FOLD_CASE);
	FastaReader sequence2(file2.c_str(), FastaReader::NO_FOLD_CASE);
	FastqRecord rec1;
	FastqRecord rec2;
	//hits results stored in hashmap of filter names and hits
	unordered_map<string, size_t> hits1(filterNum);
	unordered_map<string, size_t> hits2(filterNum);

	if (tileModifier == 0) {
		while (sequence1 >> rec1 && sequence2 >> rec2) {
			++totalReads;
			if (totalReads % 1000000 == 0) {
				cerr << "Currently Reading Read Number: " << totalReads << endl;
			}

			//initialize hits to zero
			initHits(hits1);
			initHits(hits2);

			//for each hashSigniture/kmer combo multi, cut up read into kmer sized used
			for (vector<string>::const_iterator j = hashSigs.begin();
					j != hashSigs.end(); ++j) {
				string tempStr1 = rec1.id.substr(0, rec1.id.find_last_of("/"));
				string tempStr2 = rec2.id.substr(0, rec2.id.find_last_of("/"));
				if (tempStr1 == tempStr2) {
					evaluateRead(rec1, *j, hits1);
					evaluateRead(rec2, *j, hits2);
				} else {
					cerr << "Read IDs do not match" << "\n" << tempStr1 << "\n"
							<< tempStr2 << endl;
					exit(1);
				}
			}

			string readID = rec1.id.substr(0, rec1.id.length() - 2);

			//print hit results to read status
			readStatusOutput
					<< getReadStatStrPair(readID, rec1.seq.length(),
							rec1.seq.length(), hits1, hits2);

			//Evaluate hit data and record for summary
			const string &outputFileName = resSummary.updateSummaryData(
					rec1.seq.length(), rec2.seq.length(), hits1, hits2);

			if (outputType == "fa") {
				(*outputFiles[outputFileName + "1"]) << ">" << rec1.id << "\n"
						<< rec1.seq << "\n";
				(*outputFiles[outputFileName + "2"]) << ">" << rec2.id << "\n"
						<< rec2.seq << "\n";
			} else {
				(*outputFiles[outputFileName + "1"]) << "@" << rec1.id << "\n"
						<< rec1.seq << "\n+\n" << rec1.qual << "\n";
				(*outputFiles[outputFileName + "2"]) << "@" << rec2.id << "\n"
						<< rec2.seq << "\n+\n" << rec2.qual << "\n";
			}
		}
	} else {

		while (sequence1 >> rec1 && sequence2 >> rec2) {
			++totalReads;
			if (totalReads % 1000000 == 0) {
				cerr << "Currently Reading Read Number: " << totalReads << endl;
			}

			//initialize hits to zero
			initHits(hits1);
			initHits(hits2);

			//for each hashSigniture/kmer combo multi, cut up read into kmer sized used
			for (vector<string>::const_iterator j = hashSigs.begin();
					j != hashSigs.end(); ++j) {
				string tempStr1 = rec1.id.substr(0, rec1.id.find_last_of("/"));
				string tempStr2 = rec2.id.substr(0, rec2.id.find_last_of("/"));
				if (tempStr1 == tempStr2) {
					evaluateRead(rec1, *j, hits1, tileModifier);
					evaluateRead(rec2, *j, hits2, tileModifier);
				} else {
					cerr << "Read IDs do not match" << "\n" << tempStr1 << "\n"
							<< tempStr2 << endl;
					exit(1);
				}
			}

			string readID = rec1.id.substr(0, rec1.id.length() - 2);

			//print hit results to read status
			readStatusOutput
					<< getReadStatStrPair(readID, rec1.seq.length(),
							rec1.seq.length(), hits1, hits2);

			//Evaluate hit data and record for summary
			const string &outputFileName = resSummary.updateSummaryData(
					rec1.seq.length(), rec2.seq.length(), hits1, hits2);

			if (outputType == "fa") {
				(*outputFiles[outputFileName + "1"]) << ">" << rec1.id << "\n"
						<< rec1.seq << "\n";
				(*outputFiles[outputFileName + "2"]) << ">" << rec2.id << "\n"
						<< rec2.seq << "\n";
			} else {
				(*outputFiles[outputFileName + "1"]) << "@" << rec1.id << "\n"
						<< rec1.seq << "\n+\n" << rec1.qual << "\n";
				(*outputFiles[outputFileName + "2"]) << "@" << rec2.id << "\n"
						<< rec2.seq << "\n+\n" << rec2.qual << "\n";
			}
		}
	}
	if (sequence2 >> rec2 && sequence1.eof() && sequence2.eof()) {
		cerr
				<< "error: eof bit not flipped. Input files may be different lengths"
				<< endl;
	}

	readStatusOutput.close();

	//close sorting files
	for (unordered_map<string, shared_ptr<Dynamicofstream> >::iterator j =
			outputFiles.begin(); j != outputFiles.end(); ++j) {
		j->second->close();
	}

	cerr << "Total Reads:" << totalReads << endl;

	Dynamicofstream summaryOutput(prefix + "_summary.tsv");
	summaryOutput << resSummary.getResultsSummary(totalReads);
	summaryOutput.close();

	if (maxHitValue > 0) {
		Dynamicofstream countSummaryOutput(prefix + "_rawCounts.tsv");
		countSummaryOutput << resSummary.getCountSummary(totalReads);
		countSummaryOutput.close();
	}
}

/*
 * Filters reads -> uses paired end information
 * Assumes only one hash signature exists (load only filters with same
 * hash functions)
 */
void BioBloomClassifier::filterPairBAM(const string &file) {
	Dynamicofstream readStatusOutput(prefix + "_status.tsv" + postfix);

	//results summary object
	ResultsManager resSummary(hashSigs, filters, infoFiles, minHit,
			percentMinHit, maxHitValue, tileModifier);

	unordered_map<string, FastqRecord> unPairedReads;

	size_t totalReads = 0;

	//print out header info and initialize variables for summary
	readStatusOutput << getReadSummaryHeader(hashSigs);

	cerr << "Filtering Start" << "\n";

	FastaReader sequence(file.c_str(), FastaReader::NO_FOLD_CASE);
	//hits results stored in hashmap of filter names and hits
	unordered_map<string, size_t> hits1(filterNum);
	unordered_map<string, size_t> hits2(filterNum);

	if (tileModifier == 0) {
		while (!sequence.eof()) {
			FastqRecord rec;
			if (sequence >> rec) {
				string readID = rec.id.substr(0, rec.id.length() - 2);
				if (unPairedReads.find(readID) != unPairedReads.end()) {

					const FastqRecord &rec1 =
							rec.id.at(rec.id.length() - 1) == '1' ?
									rec : unPairedReads[readID];
					const FastqRecord &rec2 =
							rec.id.at(rec.id.length() - 1) == '2' ?
									rec : unPairedReads[readID];

					++totalReads;
					if (totalReads % 1000000 == 0) {
						cerr << "Currently Reading Read Number: " << totalReads
								<< endl;
					}

					//initialize hits to zero
					initHits(hits1);
					initHits(hits2);

					//for each hashSigniture/kmer combo multi, cut up read into kmer sized used
					for (vector<string>::const_iterator j = hashSigs.begin();
							j != hashSigs.end(); ++j) {
						string tempStr1 = rec1.id.substr(0,
								rec1.id.find_last_of("/"));
						string tempStr2 = rec2.id.substr(0,
								rec2.id.find_last_of("/"));
						if (tempStr1 == tempStr2) {
							evaluateRead(rec1, *j, hits1);
							evaluateRead(rec2, *j, hits2);
						} else {
							cerr << "Read IDs do not match" << "\n" << tempStr1
									<< "\n" << tempStr2 << endl;
							exit(1);
						}
					}

					//print hit results to read status
					readStatusOutput
							<< getReadStatStrPair(readID, rec1.seq.length(),
									rec2.seq.length(), hits1, hits2);

					//Evaluate hit data and record for summary
					resSummary.updateSummaryData(rec1.seq.length(),
							rec2.seq.length(), hits1, hits2);

					//clean up reads
					unPairedReads.erase(readID);

				} else {
					unPairedReads[readID] = rec;
				}
			}
		}
	} else {
		while (!sequence.eof()) {
			FastqRecord rec;
			if (sequence >> rec) {
				string readID = rec.id.substr(0, rec.id.length() - 2);
				if (unPairedReads.find(readID) != unPairedReads.end()) {

					const FastqRecord &rec1 =
							rec.id.at(rec.id.length() - 1) == '1' ?
									rec : unPairedReads[readID];
					const FastqRecord &rec2 =
							rec.id.at(rec.id.length() - 1) == '2' ?
									rec : unPairedReads[readID];

					++totalReads;
					if (totalReads % 1000000 == 0) {
						cerr << "Currently Reading Read Number: " << totalReads
								<< endl;
					}

					//initialize hits to zero
					initHits(hits1);
					initHits(hits2);

					//for each hashSigniture/kmer combo multi, cut up read into kmer sized used
					for (vector<string>::const_iterator j = hashSigs.begin();
							j != hashSigs.end(); ++j) {
						string tempStr1 = rec1.id.substr(0,
								rec1.id.find_last_of("/"));
						string tempStr2 = rec2.id.substr(0,
								rec2.id.find_last_of("/"));
						if (tempStr1 == tempStr2) {
							evaluateRead(rec1, *j, hits1, tileModifier);
							evaluateRead(rec2, *j, hits2, tileModifier);
						} else {
							cerr << "Read IDs do not match" << "\n" << tempStr1
									<< "\n" << tempStr2 << endl;
							exit(1);
						}
					}

					//print hit results to read status
					readStatusOutput
							<< getReadStatStrPair(readID, rec1.seq.length(),
									rec2.seq.length(), hits1, hits2);

					//Evaluate hit data and record for summary
					resSummary.updateSummaryData(rec1.seq.length(),
							rec2.seq.length(), hits1, hits2);

					//clean up reads
					unPairedReads.erase(readID);

				} else {
					unPairedReads[readID] = rec;
				}
			}
		}
	}

	readStatusOutput.close();

	Dynamicofstream summaryOutput(prefix + "_summary.tsv");
	summaryOutput << resSummary.getResultsSummary(totalReads);
	summaryOutput.close();

	if (maxHitValue > 0) {
		Dynamicofstream countSummaryOutput(prefix + "_rawCounts.tsv");
		countSummaryOutput << resSummary.getCountSummary(totalReads);
		countSummaryOutput.close();
	}
}

/*
 * Filters reads -> uses paired end information
 * Assumes only one hash signature exists (load only filters with same
 * hash functions)
 * Prints reads into seperate files
 */
void BioBloomClassifier::filterPairBAMPrint(const string &file,
		const string &outputType) {
	Dynamicofstream readStatusOutput(prefix + "_status.tsv" + postfix);

	//results summary object
	ResultsManager resSummary(hashSigs, filters, infoFiles, minHit,
			percentMinHit, maxHitValue, tileModifier);

	unordered_map<string, FastqRecord> unPairedReads;

	size_t totalReads = 0;

	unordered_map<string, shared_ptr<Dynamicofstream> > outputFiles;
	shared_ptr<Dynamicofstream> noMatch1(
			new Dynamicofstream(
					prefix + "_" + noMatch + "_1." + outputType + postfix));
	shared_ptr<Dynamicofstream> noMatch2(
			new Dynamicofstream(
					prefix + "_" + noMatch + "_2." + outputType + postfix));
	shared_ptr<Dynamicofstream> multiMatch1(
			new Dynamicofstream(
					prefix + "_" + multiMatch + "_1." + outputType + postfix));
	shared_ptr<Dynamicofstream> multiMatch2(
			new Dynamicofstream(
					prefix + "_" + multiMatch + "_2." + outputType + postfix));
	outputFiles[noMatch + "1"] = noMatch1;
	outputFiles[noMatch + "2"] = noMatch2;
	outputFiles[multiMatch + "1"] = multiMatch1;
	outputFiles[multiMatch + "2"] = multiMatch2;

	//initialize variables and print filter ids
	for (vector<string>::const_iterator j = hashSigs.begin();
			j != hashSigs.end(); ++j) {
		const vector<string> idsInFilter = (*filters[*j]).getFilterIds();
		for (vector<string>::const_iterator i = idsInFilter.begin();
				i != idsInFilter.end(); ++i) {
			shared_ptr<Dynamicofstream> temp1(
					new Dynamicofstream(
							prefix + "_" + *i + "_1." + outputType + postfix));
			shared_ptr<Dynamicofstream> temp2(
					new Dynamicofstream(
							prefix + "_" + *i + "_2." + outputType + postfix));
			outputFiles[*i + "1"] = temp1;
			outputFiles[*i + "2"] = temp2;
		}
	}

	//print out header info and initialize variables for summary
	readStatusOutput << getReadSummaryHeader(hashSigs);
	cerr << "Filtering Start" << "\n";

	FastaReader sequence(file.c_str(), FastaReader::NO_FOLD_CASE);
	//hits results stored in hashmap of filter names and hits
	unordered_map<string, size_t> hits1(filterNum);
	unordered_map<string, size_t> hits2(filterNum);

	if (tileModifier == 0) {
		while (!sequence.eof()) {
			FastqRecord rec;
			if (sequence >> rec) {
				string readID = rec.id.substr(0, rec.id.length() - 2);
				if (unPairedReads.find(readID) != unPairedReads.end()) {

					const FastqRecord &rec1 =
							rec.id.at(rec.id.length() - 1) == '1' ?
									rec : unPairedReads[readID];
					const FastqRecord &rec2 =
							rec.id.at(rec.id.length() - 1) == '2' ?
									rec : unPairedReads[readID];

					++totalReads;
					if (totalReads % 1000000 == 0) {
						cerr << "Currently Reading Read Number: " << totalReads
								<< endl;
					}

					//initialize hits to zero
					initHits(hits1);
					initHits(hits2);

					//for each hashSigniture/kmer combo multi, cut up read into kmer sized used
					for (vector<string>::const_iterator j = hashSigs.begin();
							j != hashSigs.end(); ++j) {
						string tempStr1 = rec1.id.substr(0,
								rec1.id.find_last_of("/"));
						string tempStr2 = rec2.id.substr(0,
								rec2.id.find_last_of("/"));
						if (tempStr1 == tempStr2) {
							evaluateRead(rec1, *j, hits1);
							evaluateRead(rec2, *j, hits2);
						} else {
							cerr << "Read IDs do not match" << "\n" << tempStr1
									<< "\n" << tempStr2 << endl;
							exit(1);
						}
					}

					//print hit results to read status
					readStatusOutput
							<< getReadStatStrPair(readID, rec1.seq.length(),
									rec2.seq.length(), hits1, hits2);

					//Evaluate hit data and record for summary
					const string &outputFileName = resSummary.updateSummaryData(
							rec1.seq.length(), rec2.seq.length(), hits1, hits2);

					if (outputType == "fa") {
						(*outputFiles[outputFileName + "1"]) << ">" << rec1.id
								<< "\n" << rec1.seq << "\n";
						(*outputFiles[outputFileName + "2"]) << ">" << rec2.id
								<< "\n" << rec2.seq << "\n";
					} else {
						(*outputFiles[outputFileName + "1"]) << "@" << rec1.id
								<< "\n" << rec1.seq << "\n+\n" << rec1.qual
								<< "\n";
						(*outputFiles[outputFileName + "2"]) << "@" << rec2.id
								<< "\n" << rec2.seq << "\n+\n" << rec2.qual
								<< "\n";
					}

					//clean up reads
					unPairedReads.erase(readID);

				} else {
					unPairedReads[readID] = rec;
				}
			}
		}
	} else {
		while (!sequence.eof()) {
			FastqRecord rec;
			if (sequence >> rec) {
				string readID = rec.id.substr(0, rec.id.length() - 2);
				if (unPairedReads.find(readID) != unPairedReads.end()) {

					const FastqRecord &rec1 =
							rec.id.at(rec.id.length() - 1) == '1' ?
									rec : unPairedReads[readID];
					const FastqRecord &rec2 =
							rec.id.at(rec.id.length() - 1) == '2' ?
									rec : unPairedReads[readID];

					++totalReads;
					if (totalReads % 1000000 == 0) {
						cerr << "Currently Reading Read Number: " << totalReads
								<< endl;
					}

					//initialize hits to zero
					initHits(hits1);
					initHits(hits2);

					//for each hashSigniture/kmer combo multi, cut up read into kmer sized used
					for (vector<string>::const_iterator j = hashSigs.begin();
							j != hashSigs.end(); ++j) {
						string tempStr1 = rec1.id.substr(0,
								rec1.id.find_last_of("/"));
						string tempStr2 = rec2.id.substr(0,
								rec2.id.find_last_of("/"));
						if (tempStr1 == tempStr2) {
							evaluateRead(rec1, *j, hits1, tileModifier);
							evaluateRead(rec2, *j, hits2, tileModifier);
						} else {
							cerr << "Read IDs do not match" << "\n" << tempStr1
									<< "\n" << tempStr2 << endl;
							exit(1);
						}
					}

					//print hit results to read status
					readStatusOutput
							<< getReadStatStrPair(readID, rec1.seq.length(),
									rec2.seq.length(), hits1, hits2);

					//Evaluate hit data and record for summary
					const string &outputFileName = resSummary.updateSummaryData(
							rec1.seq.length(), rec2.seq.length(), hits1, hits2);

					if (outputType == "fa") {
						(*outputFiles[outputFileName + "1"]) << ">" << rec1.id
								<< "\n" << rec1.seq << "\n";
						(*outputFiles[outputFileName + "2"]) << ">" << rec2.id
								<< "\n" << rec2.seq << "\n";
					} else {
						(*outputFiles[outputFileName + "1"]) << "@" << rec1.id
								<< "\n" << rec1.seq << "\n+\n" << rec1.qual
								<< "\n";
						(*outputFiles[outputFileName + "2"]) << "@" << rec2.id
								<< "\n" << rec2.seq << "\n+\n" << rec2.qual
								<< "\n";
					}

					//clean up reads
					unPairedReads.erase(readID);

				} else {
					unPairedReads[readID] = rec;
				}
			}
		}
	}

	//close sorting files
	for (unordered_map<string, shared_ptr<Dynamicofstream> >::iterator j =
			outputFiles.begin(); j != outputFiles.end(); ++j) {
		j->second->close();
	}

	Dynamicofstream summaryOutput(prefix + "_summary.tsv");
	summaryOutput << resSummary.getResultsSummary(totalReads);
	summaryOutput.close();

	if (maxHitValue > 0) {
		Dynamicofstream countSummaryOutput(prefix + "_rawCounts.tsv");
		countSummaryOutput << resSummary.getCountSummary(totalReads);
		countSummaryOutput.close();
	}
}

//helper methods

/*
 * Loads list of filters into memory
 * todo: Implement non-block I/O when loading multiple filters at once
 */
void BioBloomClassifier::loadFilters(const vector<string> &filterFilePaths) {
	cerr << "Starting to Load Filters." << endl;
	//load up files
	for (vector<string>::const_iterator it = filterFilePaths.begin();
			it != filterFilePaths.end(); ++it) {
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
		shared_ptr<BloomFilterInfo> info(new BloomFilterInfo(infoFileName));
		//append kmer size to hash signature to insure correct kmer size is used
		stringstream hashSig;
		hashSig << info->getHashNum() << info->getKmerSize();

		//if hashSig exists add filter to list
		if (infoFiles.count(hashSig.str()) == 1) {
			infoFiles[hashSig.str()].push_back(info);
			filters[hashSig.str()]->addFilter(info->getCalcuatedFilterSize(),
					info->getFilterID(), *it);
		} else {
			vector<shared_ptr<BloomFilterInfo> > tempVect;
			tempVect.push_back(info);
			hashSigs.push_back(hashSig.str());
			shared_ptr<MultiFilter> temp(
					new MultiFilter(info->getHashNum(), info->getKmerSize()));
			filters[hashSig.str()] = temp;
			filters[hashSig.str()]->addFilter(info->getCalcuatedFilterSize(),
					info->getFilterID(), *it);
			infoFiles[hashSig.str()] = tempVect;
		}
		cerr << "Loaded Filter: " + info->getFilterID() << endl;
	}
	cerr << "Filter Loading Complete." << endl;
}

/*
 * checks if file exists
 */
const bool BioBloomClassifier::fexists(const string &filename) const {
	ifstream ifile(filename.c_str());
	return ifile;
}

/*
 * For a single read evaluate hits for a single hash signature
 * Sections with ambiguity bases are treated as misses
 * Updates hits value to number of hits (hashSig is used to as key)
 * Faster variant that assume there a redundant tile of 0
 */
void BioBloomClassifier::evaluateRead(const FastqRecord &rec,
		const string &hashSig, unordered_map<string, size_t> &hits) {
	//get filterIDs to iterate through has in a consistent order
	const vector<string> &idsInFilter = (*filters[hashSig]).getFilterIds();

	//get kmersize for set of info files
	uint16_t kmerSize = infoFiles.at(hashSig).front()->getKmerSize();

	//Establish tiling pattern
	uint16_t startModifier1 = (rec.seq.length() % kmerSize) / 2;
	size_t currentKmerNum = 0;

	ReadsProcessor proc(kmerSize);
	//cut read into kmer size given
	while (rec.seq.length() >= (currentKmerNum + 1) * kmerSize) {

		const char* currentKmer = proc.prepSeq(rec.seq,
				currentKmerNum * kmerSize + startModifier1);

		//check to see if string is invalid
		if (*currentKmer != 0) {
			const unordered_map<string, bool> &results =
					filters[hashSig]->multiContains(currentKmer);

			//record hit number in order
			for (vector<string>::const_iterator i = idsInFilter.begin();
					i != idsInFilter.end(); ++i) {
				if (results.find(*i)->second) {
					++hits[*i];
				}
			}
		}
		++currentKmerNum;
	}
}

/*
 * Assumes single filter is used
 */
size_t BioBloomClassifier::evaluateReadSingle(const FastqRecord &rec,
		const BloomFilter &filter) {
	//break read into tiles
	//evaluate read

	//depending on results evaluate pairs 3 extra tiles per pair, at median points
	//if only single section passes, evaluate read 4 times, shift only one tile
}

/*
 * For a single read evaluate hits for a single hash signature
 * Sections with ambiguity bases are treated as misses
 * Updates hits value to number of hits (hashSig is used to as key)
 */
void BioBloomClassifier::evaluateRead(const FastqRecord &rec,
		const string &hashSig, unordered_map<string, size_t> &hits,
		uint8_t redundantRead) {
	//get filterIDs to iterate through has in a consistent order
	const vector<string> &idsInFilter = (*filters[hashSig]).getFilterIds();

	//get kmersize for set of info files
	uint16_t kmerSize = infoFiles.at(hashSig).front()->getKmerSize();

	//Establish tiling pattern
	uint16_t startModifier = (rec.seq.length() % (kmerSize)) / 2;

	ReadsProcessor proc(kmerSize);

	size_t currentLoc = 0;
	//cut read into kmer size + tilemodifier given
	while (rec.seq.length() >= currentLoc + kmerSize) {

		unordered_map<string, bool> tempResults;
		vector<string> filterIDs = idsInFilter;

		for (vector<string>::const_iterator i = idsInFilter.begin();
				i != idsInFilter.end(); ++i) {
			tempResults[*i] = true;
		}

		//shift right
		if (rec.seq.length() != currentLoc + kmerSize) {
			for (uint8_t j = 0; j < redundantRead; ++j) {

				if (filterIDs.size() == 0) {
					break;
				}

				const char* currentKmer = proc.prepSeq(rec.seq,
						currentLoc + startModifier + j);
				vector<string> tempFilterIDs;

				//check to see if string is invalid
				if (*currentKmer != 0) {

					const unordered_map<string, bool> &results =
							filters[hashSig]->multiContains(currentKmer,
									filterIDs);

					for (vector<string>::iterator i = filterIDs.begin();
							i != filterIDs.end(); ++i) {
						if (!results.find(*i)->second) {
							tempResults[*i] = false;
						} else {
							tempFilterIDs.push_back(*i);
						}
					}
				} else {
					for (vector<string>::iterator i = filterIDs.begin();
							i != filterIDs.end(); ++i) {
						tempResults[*i] = false;
					}
					break;
				}
				filterIDs = tempFilterIDs;
			}
		}
		//shift left on last read segment
		else {
			for (uint8_t j = 0; j < redundantRead; ++j) {

				if (filterIDs.size() == 0) {
					break;
				}

				const char* currentKmer = proc.prepSeq(rec.seq,
						currentLoc + startModifier - 1);
				vector<string> tempFilterIDs;

				//check to see if string is invalid
				if (*currentKmer != 0) {

					const unordered_map<string, bool> &results =
							filters[hashSig]->multiContains(currentKmer,
									filterIDs);

					for (vector<string>::iterator i = filterIDs.begin();
							i != filterIDs.end(); ++i) {
						if (!results.find(*i)->second) {
							tempResults[*i] = false;
						} else {
							tempFilterIDs.push_back(*i);
						}
					}
				} else {
					for (vector<string>::iterator i = filterIDs.begin();
							i != filterIDs.end(); ++i) {
						tempResults[*i] = false;
					}
					break;
				}
				filterIDs = tempFilterIDs;
			}
		}

		//record hit number in order
		for (vector<string>::const_iterator i = idsInFilter.begin();
				i != idsInFilter.end(); ++i) {
			if (tempResults.find(*i)->second) {
				++hits[*i];
			}
		}
		currentLoc += kmerSize;
	}
}

///*
// * For a single read evaluate hits for a single hash signature
// * Sections with ambiguity bases are treated as misses
// * Updates hits value to number of hits (hashSig is used to as key)
// */
//void BioBloomClassifier::evaluateRead(const FastqRecord &rec,
//		const string &hashSig, unordered_map<string, size_t> &hits,
//		uint8_t redundantRead)
//{
//	//get filterIDs to iterate through has in a consistent order
//	const vector<string> &idsInFilter = (*filters[hashSig]).getFilterIds();
//
//	//get kmersize for set of info files
//	uint16_t kmerSize = infoFiles.at(hashSig).front()->getKmerSize();
//
//	//Establish tiling pattern
//	uint16_t startModifier = (rec.seq.length() % (kmerSize + redundantRead)) / 2;
//
//	ReadsProcessor proc(kmerSize);
//
//	size_t currentLoc = 0;
//	//cut read into kmer size + tilemodifier given
//	while (rec.seq.length() >= currentLoc + kmerSize + redundantRead) {
//
//		unordered_map<string, bool> tempResults;
//		vector<string> filterIDs = idsInFilter;
//
//		for (vector<string>::const_iterator i = idsInFilter.begin();
//				i != idsInFilter.end(); ++i)
//		{
//			tempResults[*i] = true;
//		}
//
//		for (uint8_t j = 0; j <= redundantRead; ++j) {
//
//			const string &currentKmer = proc.prepSeq(rec.seq,
//					currentLoc + startModifier + j);
//			vector<string> tempFilterIDs;
//
//			//check to see if string is invalid
//			if (!currentKmer.empty()) {
//
//				const unordered_map<string, bool> &results =
//						filters[hashSig]->multiContains(currentKmer, filterIDs);
//
//				for (vector<string>::iterator i = filterIDs.begin();
//						i != filterIDs.end(); ++i)
//				{
//					if (!results.find(*i)->second) {
//						tempResults[*i] = false;
//					} else {
//						tempFilterIDs.push_back(*i);
//					}
//				}
//			} else {
//				for (vector<string>::iterator i = filterIDs.begin();
//						i != filterIDs.end(); ++i)
//				{
//					tempResults[*i] = false;
//				}
//			}
//			filterIDs = tempFilterIDs;
//		}
//
//		//record hit number in order
//		for (vector<string>::const_iterator i = idsInFilter.begin();
//				i != idsInFilter.end(); ++i)
//		{
//			if (tempResults.find(*i)->second) {
//				++hits[*i];
//			}
//		}
//		currentLoc += kmerSize + redundantRead;
//	}
//}

/*
 * Initializes Summary Variables. Also prints heads for read status.
 */
const string BioBloomClassifier::getReadSummaryHeader(
		const vector<string> &hashSigs) {
	stringstream readStatusOutput;
	readStatusOutput << "readID\tseqSize";

	//initialize variables and print filter ids
	for (vector<string>::const_iterator j = hashSigs.begin();
			j != hashSigs.end(); ++j) {
		vector<string> idsInFilter = (*filters[*j]).getFilterIds();
		for (vector<string>::const_iterator i = idsInFilter.begin();
				i != idsInFilter.end(); ++i) {
			readStatusOutput << "\t" << *i << "_"
					<< (*(infoFiles[*j].front())).getKmerSize() << "_"
					<< tileModifier;
		}
	}
	readStatusOutput << "\n";
	return readStatusOutput.str();
}

/*
 * Initializes hits results to zero
 */
void BioBloomClassifier::initHits(unordered_map<string, size_t> &hits) {
	//initialize hits to zero
	for (vector<string>::const_iterator j = hashSigs.begin();
			j != hashSigs.end(); ++j) {
		const vector<string> &idsInFilter = (*filters[*j]).getFilterIds();
		for (vector<string>::const_iterator i = idsInFilter.begin();
				i != idsInFilter.end(); ++i) {
			hits[*i] = 0;
		}
	}
}

/*
 * return results of hits to be output for read status
 */
const string BioBloomClassifier::getReadStatStr(string const &readID,
		size_t readLength, unordered_map<string, size_t> &hits) {
	stringstream str;
	str << readID << "\t" << readLength;
	//print readID
	for (vector<string>::const_iterator j = hashSigs.begin();
			j != hashSigs.end(); ++j) {
		//update summary
		const vector<string> &idsInFilter = (*filters[*j]).getFilterIds();
		for (vector<string>::const_iterator i = idsInFilter.begin();
				i != idsInFilter.end(); ++i) {
			//print to file
			str << "\t" << hits[*i];
		}
	}
	str << "\n";
	return str.str();
}

/*
 * return results of hits to be output for read status for paired end mode
 */
const string BioBloomClassifier::getReadStatStrPair(string const &readID,
		size_t readLength1, size_t readLength2,
		unordered_map<string, size_t> &hits1,
		unordered_map<string, size_t> &hits2) {
	stringstream str;
	str << readID << "\t" << readLength1 << "|" << readLength2;
	//print readID
	for (vector<string>::const_iterator j = hashSigs.begin();
			j != hashSigs.end(); ++j) {
		//update summary
		const vector<string> &idsInFilter = (*filters[*j]).getFilterIds();
		for (vector<string>::const_iterator i = idsInFilter.begin();
				i != idsInFilter.end(); ++i) {
			//print to file
			str << "\t" << hits1[*i] << "|" << hits2[*i];
		}
	}
	str << "\n";
	return str.str();
}

BioBloomClassifier::~BioBloomClassifier() {
}

