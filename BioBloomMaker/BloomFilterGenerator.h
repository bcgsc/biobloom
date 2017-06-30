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
#include "Common/BloomFilter.h"
#include "Common/SeqEval.h"
#include "Common/SeqEval.h"
#include "DataLayer/kseq.h"
#include <zlib.h>
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
			const SeqEval::EvalMode evalMode, bool printReads);
	size_t generateProgressive(const string &filename, double score,
			const string &file1, const string &file2, createMode mode,
			const SeqEval::EvalMode evalMode, bool printReads,
			const string &subtractFilter);
	size_t generateProgressiveBait(const string &filename, double score,
			const string &file1, const string &file2, createMode mode,
			const SeqEval::EvalMode evalMode, bool printReads,
			const string &subtractFilter);

	size_t generateProgressive(const string &filename, double score,
			const vector<string> &files1, const vector<string> &files2,
			createMode mode, const SeqEval::EvalMode evalMode, bool printReads,
			const string &subtractFilter);
	size_t generateProgressiveBait(const string &filename, double score,
			const vector<string> &files1, const vector<string> &files2,
			createMode mode, const SeqEval::EvalMode evalMode, bool printReads,
			const string &subtractFilter);

	size_t generateProgressive(const string &filename, double score,
			const vector<string> &files, const SeqEval::EvalMode evalMode,
			bool printReads, const string &subtractFilter);

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

	inline size_t loadFilterFast(BloomFilter &bf, size_t &totalEntries) {
		unsigned threadNum = omp_get_max_threads();
		vector<boost::shared_ptr<ReadsProcessor> > procs(threadNum);
		vector<boost::shared_ptr<vector<size_t> > > tempHashValues(threadNum);

		size_t redundancy = 0;

		//each thread gets its own thread processor
		for (unsigned i = 0; i < threadNum; ++i) {
			procs[i] = boost::shared_ptr<ReadsProcessor>(
					new ReadsProcessor(m_kmerSize));
			tempHashValues[i] = boost::shared_ptr<vector<size_t> >(
					new vector<size_t>(m_hashNum));
		}

		int kmerSize = m_kmerSize;

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
					int screeningLoc = 0;
					//k-merize and insert into bloom filter
					while ((screeningLoc + kmerSize) <= l) {
						const unsigned char* currentKmer =
								procs[omp_get_thread_num()]->prepSeq(tempStr,
										screeningLoc);
						if (currentKmer != NULL) {
							bool found =
									bf.insertAndCheck(
											multiHash(currentKmer, m_hashNum,
													m_kmerSize,
													*tempHashValues[omp_get_thread_num()]));
							tempRedund += found;
							tempTotal += !found;
						}
						++screeningLoc;
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

	inline size_t loadFilterFastSubtract(BloomFilter &bf, BloomFilter &bfsub,
			size_t &totalEntries) {
		unsigned threadNum = omp_get_max_threads();
		vector<boost::shared_ptr<ReadsProcessor> > procs(threadNum);
		vector<boost::shared_ptr<vector<size_t> > > tempHashValues(threadNum);

		size_t kmerRemoved = 0;
		size_t redundancy = 0;

		//each thread gets its own thread processor
		for (unsigned i = 0; i < threadNum; ++i) {
			procs[i] = boost::shared_ptr<ReadsProcessor>(
					new ReadsProcessor(m_kmerSize));
			tempHashValues[i] = boost::shared_ptr<vector<size_t> >(
					new vector<size_t>(m_hashNum));
		}

		int kmerSize = m_kmerSize;

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
					int screeningLoc = 0;
					//k-merize and insert into bloom filter
					while ((screeningLoc + kmerSize) <= l) {
						const unsigned char* currentKmer =
								procs[omp_get_thread_num()]->prepSeq(tempStr,
										screeningLoc);
						if (currentKmer != NULL) {
							if (bfsub.contains(currentKmer)) {
								++kmerRemoved;
							} else {
								bool found =
										bf.insertAndCheck(
												multiHash(currentKmer,
														m_hashNum, m_kmerSize,
														*tempHashValues[omp_get_thread_num()]));
								tempRedund += found;
								tempTotal += !found;
							}
						}
						++screeningLoc;
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

	inline void insertKmer(const vector<size_t> &currentSeq,
			BloomFilter &filter) {
#pragma omp atomic
		m_totalEntries += !filter.insertAndCheck(currentSeq);
	}

	inline void insertKmer(const unsigned char* currentSeq,
			BloomFilter &filter) {
		if (currentSeq != NULL) {
#pragma omp atomic
			m_totalEntries += !filter.insertAndCheck(
					multiHash(currentSeq, m_hashNum, m_kmerSize));
		}
	}
};

#endif /* BLOOMFILTERGENERATOR_H_ */
