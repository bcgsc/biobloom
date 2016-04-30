/*
 * BloomMapGenerator.cpp
 *
 *  Created on: Mar 17, 2016
 *      Author: cjustin
 */

#include "BloomMapGenerator.h"
#include <zlib.h>
#include <stdio.h>
#include "DataLayer/kseq.h"
KSEQ_INIT(gzFile, gzread)

/*
 * If size is set to 0 (or not set) a filter size will be estimated based
 * on the number of k-mers in the file
 */
BloomMapGenerator::BloomMapGenerator(vector<string> const &filenames,
		unsigned kmerSize, size_t numElements = 0) :
		m_kmerSize(kmerSize), m_expectedEntries(numElements), m_totalEntries(0), m_fileNames(
				filenames){
	//Instantiate dense hash map
	m_headerIDs.set_empty_key("");

	ID value = 0;
	size_t counts = 0;

	if (opt::idByFile) {
		for (unsigned i = 0; i < m_fileNames.size(); ++i) {
			gzFile fp;
			fp = gzopen(m_fileNames[i].c_str(), "r");
			kseq_t *seq = kseq_init(fp);
			m_headerIDs[m_fileNames[i].substr(m_fileNames[i].find_last_of("/")+1)] = ++value;
			int l;
			for (;;) {
				l = kseq_read(seq);
				if (l >= 0) {
					counts += seq->seq.l - m_kmerSize + 1;
				} else {
					kseq_destroy(seq);
					break;
				}
			}
			gzclose(fp);
		}
	} else {
		for (unsigned i = 0; i < m_fileNames.size(); ++i) {
			gzFile fp;
			fp = gzopen(m_fileNames[i].c_str(), "r");
			kseq_t *seq = kseq_init(fp);
			int l;
			for (;;) {
				l = kseq_read(seq);
				if (l >= 0) {
					m_headerIDs[seq->name.s] = ++value;
					counts += seq->seq.l - m_kmerSize + 1;
				} else {
					kseq_destroy(seq);
					break;
				}
			}
			gzclose(fp);
		}
	}

	//estimate number of k-mers
	if (m_expectedEntries == 0) {
		m_expectedEntries = counts;
	}
	//assume each file one line per file
	cerr << "Expected number of elements: " << m_expectedEntries << endl;
}

/* Generate the bloom filter to the output filename
 */
void BloomMapGenerator::generate(const string &filePrefix, double fpr) {
	vector<vector<unsigned> > ssVal = parseSeedString(opt::sseeds);
	size_t filterSize = calcOptimalSize(m_expectedEntries, opt::sseeds.size(), fpr);

	cerr << "bitVector Size: " << filterSize << endl;

	sdsl::bit_vector *bv = new sdsl::bit_vector(filterSize);

	cerr << "Populating initial bit vector" << endl;

	size_t uniqueCounts = 0;

	//populate sdsl bitvector (bloomFilter)
#pragma omp parallel for
	for (unsigned i = 0; i < m_fileNames.size(); ++i) {
		gzFile fp;
		fp = gzopen(m_fileNames[i].c_str(), "r");
		kseq_t *seq = kseq_init(fp);
		int l;
		for (;;) {
			l = kseq_read(seq);
			if (l >= 0) {
				//k-merize with rolling hash insert into bloom map
#pragma omp atomic
				uniqueCounts += seq->seq.l - m_kmerSize + 1 - loadSeq(*bv, string(seq->seq.s, seq->seq.l), ssVal);
			} else {
				kseq_destroy(seq);
				break;
			}
		}
		gzclose(fp);
	}

	cerr << "Approximate number of unique entries in filter: " << uniqueCounts << endl;

	cerr << "Converting bit vector to rank interleaved form" << endl;

	//init bloom map
	BloomMapSSBitVec<ID> bloomMapBV(m_expectedEntries, fpr, opt::sseeds, *bv);
	//free memory of old bv
	delete(bv);

	cerr << "Populating values of bloom map" << endl;

	//populate values in bitvector (bloomFilter)
	if (opt::idByFile) {
#pragma omp parallel for
		for (unsigned i = 0; i < m_fileNames.size(); ++i) {
			gzFile fp;
			fp = gzopen(m_fileNames[i].c_str(), "r");
			kseq_t *seq = kseq_init(fp);
			int l;
			for (;;) {
				l = kseq_read(seq);
				if (l >= 0) {
					//k-merize with rolling hash insert into bloom map
					loadSeq(bloomMapBV, string(seq->seq.s, seq->seq.l),
							m_headerIDs[m_fileNames[i].substr(m_fileNames[i].find_last_of("/")+1)]);
				} else {
					kseq_destroy(seq);
					break;
				}
			}
			gzclose(fp);
		}
	} else {
#pragma omp parallel for
		for (unsigned i = 0; i < m_fileNames.size(); ++i) {
			gzFile fp;
			fp = gzopen(m_fileNames[i].c_str(), "r");
			kseq_t *seq = kseq_init(fp);
			int l;
			for (;;) {
				l = kseq_read(seq);
				if (l >= 0) {
					//k-merize with rolling hash insert into bloom map
					loadSeq(bloomMapBV, string(seq->seq.s, seq->seq.l),
							m_headerIDs[seq->name.s]);
				} else {
					kseq_destroy(seq);
					break;
				}
			}
			gzclose(fp);
		}
	}

	bloomMapBV.setUnique(uniqueCounts);

	cerr << "Storing filter" << endl;

	//save filter
	bloomMapBV.store(filePrefix + ".bf");
	writeIDs(filePrefix + "_ids.txt", m_headerIDs);
}



BloomMapGenerator::~BloomMapGenerator() {
}

