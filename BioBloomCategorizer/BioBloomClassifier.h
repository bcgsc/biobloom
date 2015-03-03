/*
 * BioBloomClassifier.h
 *
 *  Created on: Oct 17, 2012
 *      Author: cjustin
 */

#ifndef BIOBLOOMCLASSIFIER_H_
#define BIOBLOOMCLASSIFIER_H_
#include <vector>
#include <string>
#include "boost/unordered/unordered_map.hpp"
#include "boost/shared_ptr.hpp"
#include "Common/BloomFilterInfo.h"
#include "MultiFilter.h"
#include "DataLayer/FastaReader.h"
#include "Common/ReadsProcessor.h"
#include "Common/Uncompress.h"
#include "Common/BloomFilter.h"
#include "ResultsManager.h"
#include "Common/Dynamicofstream.h"

using namespace std;
using namespace boost;

static const string NO_MATCH = "noMatch";
static const string MULTI_MATCH = "multiMatch";

/** for modes of filtering */
enum mode { COLLAB, MINHITONLY, BESTHIT, STD };

//TODO: some inlining may help performance

class BioBloomClassifier {
public:
	explicit BioBloomClassifier(const vector<string> &filterFilePaths,
			double scoreThreshold, const string &outputPrefix,
			const string &outputPostFix, unsigned streakThreshold,
			unsigned minHit, bool minHitOnly);
	void filter(const vector<string> &inputFiles);
	void filterPrint(const vector<string> &inputFiles,
			const string &outputType);
	void filterPair(const string &file1, const string &file2);
	void filterPairPrint(const string &file1, const string &file2,
			const string &outputType);
	void filterPairBAM(const string &file);
	void filterPairBAMPrint(const string &file, const string &outputType);

	void setCollabFilter()
	{
		m_mode = COLLAB;
		if (m_hashSigs.size() != 1) {
			cerr
					<< "To use collaborative filtering all filters must use the same k and same number of hash functions."
					<< endl;
			exit(1);
		}
	}

	void setInclusive()
	{
		m_inclusive = true;
	}

	void setMainFilter(const string &filtername);

	virtual ~BioBloomClassifier();

private:
	//group filters with same hash number
	unordered_map<string, vector<boost::shared_ptr<BloomFilterInfo> > > m_infoFiles;
	unordered_map<string, boost::shared_ptr<MultiFilter> > m_filters;
	unordered_map<string, boost::shared_ptr<BloomFilter> > m_filtersSingle;
	vector<string> m_filterOrder;
	vector<string> m_hashSigs;
	double m_scoreThreshold;
	unsigned m_filterNum;
	const string &m_prefix;
	const string &m_postfix;
	const unsigned m_streakThreshold;
	const unsigned m_minHit;

	// modes of filtering
	mode m_mode;

	string m_mainFilter;
	bool m_inclusive;

	void loadFilters(const vector<string> &filterFilePaths);
	bool fexists(const string &filename) const;
	void evaluateReadStd(const FastqRecord &rec, const string &hashSig,
			unordered_map<string, bool> &hits);
	void evaluateReadMin(const FastqRecord &rec, const string &hashSig,
			unordered_map<string, bool> &hits);
	void evaluateReadCollab(const FastqRecord &rec, const string &hashSig,
			unordered_map<string, bool> &hits);
	double evaluateReadBestHit(const FastqRecord &rec, const string &hashSig,
			unordered_map<string, bool> &hits);

	inline void printSingle(const FastqRecord &rec,
			unordered_map<string, bool> &hits, double score)
	{
		if (m_mainFilter != "" && hits.at(m_mainFilter)) {
			if (m_mode == BESTHIT) {
#pragma omp critical(cout)
				{
					cout << "@" << rec.id << " " << score << "\n" << rec.seq
							<< "\n+\n" << rec.qual << "\n";
				}
			} else {
#pragma omp critical(cout)
				{
					cout << rec;
				}
			}
		}
	}

	inline void printSingleToFile(const string &outputFileName,
			const FastqRecord &rec,
			unordered_map<string, boost::shared_ptr<Dynamicofstream> > &outputFiles,
			const string &outputType, double score)
	{
		if (outputType == "fa") {
			if (m_mode == BESTHIT) {
#pragma omp critical(outputFiles)
				{
					(*outputFiles[outputFileName]) << ">" << rec.id << " "
							<< score << "\n" << rec.seq << "\n";
				}
			} else {
#pragma omp critical(outputFiles)
				{
					(*outputFiles[outputFileName]) << ">" << rec.id << "\n"
							<< rec.seq << "\n";
				}
			}
		} else {
			if (m_mode == BESTHIT) {
#pragma omp critical(outputFiles)
				{
					(*outputFiles[outputFileName]) << "@" << rec.id << " "
												<< score << "\n" << rec.seq << "\n+\n" << rec.qual
												<< "\n";
				}
			} else {
#pragma omp critical(outputFiles)
				{
					(*outputFiles[outputFileName]) << "@" << rec.id << "\n"
							<< rec.seq << "\n+\n" << rec.qual << "\n";
				}
			}
		}
	}

	inline void printPair(const FastqRecord &rec1, const FastqRecord &rec2,
			unordered_map<string, bool> &hits1,
			unordered_map<string, bool> &hits2)
	{
		if (m_inclusive) {
			if (m_mainFilter != ""
					&& (hits1.at(m_mainFilter) || hits2.at(m_mainFilter)))
			{
				cout << rec1;
				cout << rec2;
			}
		} else {
			if (m_mainFilter != "" && hits1.at(m_mainFilter)
					&& hits2.at(m_mainFilter))
			{
				cout << rec1;
				cout << rec2;
			}
		}
	}

	inline void printPairToFile(const string &outputFileName,
			const FastqRecord &rec1, const FastqRecord &rec2,
			unordered_map<string, boost::shared_ptr<Dynamicofstream> > &outputFiles,
			const string &outputType)
	{
		if (outputType == "fa") {
#pragma omp critical(outputFiles)
			{
				(*outputFiles[outputFileName + "1"]) << ">" << rec1.id << "\n"
						<< rec1.seq << "\n";
				(*outputFiles[outputFileName + "2"]) << ">" << rec2.id << "\n"
						<< rec2.seq << "\n";
			}
		} else {
#pragma omp critical(outputFiles)
			{
				(*outputFiles[outputFileName + "1"]) << "@" << rec1.id << "\n"
						<< rec1.seq << "\n+\n" << rec1.qual << "\n";
				(*outputFiles[outputFileName + "2"]) << "@" << rec2.id << "\n"
						<< rec2.seq << "\n+\n" << rec2.qual << "\n";
			}
		}
	}

	inline void evaluateRead(const FastqRecord &rec, const string &hashSig,
			unordered_map<string, bool> &hits, double &score)
	{
		switch(m_mode) {
		case COLLAB:{
			evaluateReadCollab(rec, hashSig, hits);
			break;
		}
		case MINHITONLY: {
			evaluateReadMin(rec, hashSig, hits);
			break;
		}
		case BESTHIT: {
			score = evaluateReadBestHit(rec, hashSig, hits);
			break;
		}
		default: {
			evaluateReadStd(rec, hashSig, hits);
			break;
		}
		}
	}
};

#endif /* BIOBLOOMCLASSIFIER_H_ */
