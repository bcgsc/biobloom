/*
 * BloomFilterGenerator.h
 *
 *  Created on: Jun 20, 2012
 *      Author: cjustin
 */

#ifndef BLOOMFILTERGENERATOR_H_
#define BLOOMFILTERGENERATOR_H_
#include <boost/shared_ptr.hpp>
#include <vector>
#include "btl_bloomfilter/BloomFilter.hpp"
#include "btl_bloomfilter/ntHashIterator.hpp"
#include "Common/SeqEval.h"
#include "DataLayer/kseq.h"
#include <iostream>
#include <zlib.h>
#include <omp.h>
#include "Common/Options.h"
KSEQ_INIT(gzFile, gzread)

using namespace std;

enum createMode {
	PROG_STD, PROG_INC
};

class BloomFilterGenerator {
public:
	explicit BloomFilterGenerator(vector<string> const &filenames,
			unsigned kmerSize, unsigned hashNum, size_t numElements);

	explicit BloomFilterGenerator(vector<string> const &filenames,
			unsigned kmerSize, unsigned hashNum);

	size_t generate(const string &filename);
	size_t generate(const string &filename, const string &subtractFilter);
	size_t generateProgressive(const string &filename, double score,
			const string &file1, const string &file2, createMode mode,
			bool printReads);
	size_t generateProgressive(const string &filename, double score,
			const string &file1, const string &file2, createMode mode,
			bool printReads, const string &subtractFilter);
	size_t generateProgressiveBait(const string &filename, double score,
			const string &file1, const string &file2, createMode mode,
			bool printReads, const string &subtractFilter);

	size_t generateProgressive(const string &filename, double score,
			const vector<string> &files1, const vector<string> &files2,
			createMode mode, bool printReads, const string &subtractFilter);
	size_t generateProgressiveBait(const string &filename, double score,
			const vector<string> &files1, const vector<string> &files2,
			createMode mode, bool printReads, const string &subtractFilter);

	size_t generateProgressive(const string &filename, double score,
			const vector<string> &files, bool printReads,
			const string &subtractFilter);

	void setFilterSize(size_t bits);

	inline void printReadPair(const string &rec1, const string &header1,
			const string &rec2, const string &header2);
	inline void printRead(const string &rec, const string &header);
	void setHashFuncs(unsigned numFunc);
	size_t getTotalEntries() const;
	size_t getExpectedEntries() const;

	virtual ~BloomFilterGenerator();
private:
	vector<string> m_fileNames;
	unsigned m_kmerSize;
	unsigned m_hashNum;
	size_t m_expectedEntries;
	size_t m_filterSize;
	size_t m_totalEntries;

	inline size_t calcExpectedEntries() {
		size_t expectedEntries = 0;
		for (unsigned i = 0; i < m_fileNames.size(); ++i) {
			gzFile fp;
			fp = gzopen(m_fileNames[i].c_str(), "r");
			if (fp == Z_NULL) {
				cerr << "file " << m_fileNames[i] << " cannot be opened"
						<< endl;
				exit(1);
			} else {
				cerr << "Opening File " << m_fileNames[i] << endl;
			}
			kseq_t *seq = kseq_init(fp);
			int l;
			size_t length;
#pragma omp parallel private(l, length)
			for (;;) {
#pragma omp critical(kseq_read)
				{
					l = kseq_read(seq);
					length = seq->seq.l;
				}
				if (l >= 0) {
#pragma omp atomic
					expectedEntries += length - m_kmerSize + 1;
				} else {
					break;
				}
			}
			kseq_destroy(seq);
			gzclose(fp);
		}
		return (expectedEntries);
	}

	inline size_t loadFilter(BloomFilter &bf, size_t &totalEntries) {
		size_t redundancy = 0;
		for (unsigned i = 0; i < m_fileNames.size(); ++i) {
			gzFile fp;
			fp = gzopen(m_fileNames[i].c_str(), "r");
			if (fp == Z_NULL) {
				cerr << "file " << m_fileNames[i] << " cannot be opened"
						<< endl;
				exit(1);
			}
			kseq_t *seq = kseq_init(fp);
			int l;
			char * tempStr;
#pragma omp parallel private(l, tempStr)
			for (;;) {
#pragma omp critical(kseq_read)
				{
					l = kseq_read(seq);
					if (l >= 0) {
						tempStr = new char[l + 1];
						strcpy(tempStr, seq->seq.s);
					}
				}
				size_t tempRedund = 0;
				size_t tempTotal = 0;
				if (l >= 0) {
					//k-merize and insert into bloom filter
					for (ntHashIterator itr(tempStr, m_hashNum, m_kmerSize); itr != itr.end(); ++itr) {
						bool found = bf.insertAndCheck(*itr);
						tempRedund += found;
						tempTotal += !found;
					}
#pragma omp atomic
					redundancy += tempRedund;
#pragma omp atomic
					totalEntries += tempTotal;
					delete[] tempStr;
				} else {
					break;
				}
			}
			kseq_destroy(seq);
			gzclose(fp);
		}
		return redundancy;
	}

	inline void loadFilter(BloomFilter &bf, const string &str) {
		size_t tempTotal = 0;
		for (ntHashIterator itr(str, m_hashNum, m_kmerSize); itr != itr.end(); ++itr) {
			tempTotal += !bf.insertAndCheck(*itr);
		}
#pragma omp atomic
		m_totalEntries += tempTotal;
	}

	inline size_t loadFilterSubtract(BloomFilter &bf, BloomFilter &bfsub,
			size_t &totalEntries) {
		size_t kmerRemoved = 0;
		size_t redundancy = 0;

		if (bf.getHashNum() != bfsub.getHashNum()) {
			cerr << "Error: Subtraction filter's hash number "
					<< bf.getHashNum()
					<< " is a different size from than output filter's hash number "
					<< bfsub.getHashNum() << endl;
			exit(1);
		}

		if (bf.getKmerSize() != bfsub.getKmerSize()) {
			cerr << "Error: Subtraction filter's k-mer size "
					<< bf.getKmerSize()
					<< " is a different size from than output filter's k-mer size "
					<< bfsub.getKmerSize() << endl;
			exit(1);
		}

		for (unsigned i = 0; i < m_fileNames.size(); ++i) {
			gzFile fp;
			fp = gzopen(m_fileNames[i].c_str(), "r");
			if (fp == NULL) {
				cerr << "file " << m_fileNames[i] << " cannot be opened"
						<< endl;
				exit(1);
			}
			kseq_t *seq = kseq_init(fp);
			int l;
			char * tempStr;
#pragma omp parallel private(l, tempStr)
			for (;;) {
#pragma omp critical(kseq_read)
				{
					l = kseq_read(seq);
					if (l >= 0) {
						tempStr = new char[l + 1];
						strcpy(tempStr, seq->seq.s);
					}
				}
				size_t tempRedund = 0;
				size_t tempTotal = 0;
				if (l >= 0) {
					for (ntHashIterator itr(tempStr, m_hashNum, m_kmerSize); itr != itr.end(); ++itr) {
						if (bfsub.contains(*itr)) {
							++kmerRemoved;
						} else {
							bool found = bf.insertAndCheck(*itr);
							tempRedund += found;
							tempTotal += !found;
						}
					}
#pragma omp atomic
					redundancy += tempRedund;
#pragma omp atomic
					totalEntries += tempTotal;
					delete[] tempStr;
				} else {
					break;
				}
			}
			kseq_destroy(seq);
			gzclose(fp);
		}
		cerr << "Total Number of K-mers not added: " << kmerRemoved << endl;
		return redundancy;
	}

	inline void insertKmer(const size_t precomputed[], BloomFilter &filter) {
#pragma omp atomic
		m_totalEntries += !filter.insertAndCheck(precomputed);
	}

	inline void insertKmer(const size_t precomputed[], BloomFilter &filter,
			const BloomFilter &filterSub) {
		if (!opt::noRep || !filterSub.contains(precomputed)) {
#pragma omp atomic
			m_totalEntries += !filter.insertAndCheck(precomputed);
		}
	}
};

#endif /* BLOOMFILTERGENERATOR_H_ */
