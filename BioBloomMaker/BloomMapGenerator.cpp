/*
 * BloomMapGenerator.cpp
 *
 *  Created on: Mar 17, 2016
 *      Author: cjustin
 */

#include "BloomMapGenerator.h"
#include "../DataLayer/FastaReader.h"
#include "../bloomfilter/BloomMap.hpp"
#include "../bloomfilter/RollingHashIterator.h"

typedef uint16_t ID;

/*
 * If size is set to 0 (or not set) a filter size will be estimated based
 * on the number of k-mers in the file
 */
BloomMapGenerator::BloomMapGenerator(vector<string> const &filenames,
		unsigned kmerSize, unsigned hashNum, size_t numElements = 0) :
		m_kmerSize(kmerSize), m_hashNum(hashNum),  m_expectedEntries(
				numElements), m_totalEntries(0), m_fileNames(filenames) {
	//estimate number of k-mers
	if (numElements == 0) {
		//assume each file one line per file
#pragma omp parallel
		for (vector<string>::const_iterator it = m_fileNames.begin();
				it != m_fileNames.end(); ++it) {
			FastaReader sequence(it->c_str(), FastaReader::NO_FOLD_CASE);
			for (FastqRecord rec;;) {
				bool good;
#pragma omp critical(sequence)
				{
					good = sequence >> rec;
				}
				if (good) {
					m_expectedEntries += rec.seq.length() - m_kmerSize;
				}
			}
		}
	}
}

void loadSeq(BloomMap<ID> bloomMap, unsigned k, unsigned hashNum, const string& seq, ID value)
	{
		if (seq.size() < k)
			return;
		/* init rolling hash state and compute hash values for first k-mer */
		RollingHashIterator itr(seq, hashNum, k);
		while (itr != itr.end()) {
			bloomMap.insert(*itr, value);
			itr++;
		}
	}

/* Generate the bloom filter to the input filename
 * Returns number of redundant entries
 */
size_t BloomMapGenerator::generate(const string &filename) {
	//init bloom map
	//for each file
	//read file, k-merize with rolling hash insert into bloom map
	//assign header to unique ID
	//if collision add new entry
	//save filter
	
	bool verbose = true;
	uint64_t count = 0;
	BloomMap<ID> bloomMap(m_expectedEntries, m_hashNum, m_kmerSize);
	ID value = 0;
 
	assert(!filename.empty());
	if (verbose)
		std::cerr << "Reading `" << filename << "'...\n";
	FastaReader in(filename.c_str(), FastaReader::NO_FOLD_CASE); //changed from FOLD_CASE
	for (FastqRecord rec;;) {
		bool good;
		good = in >> rec;
		if (good) {
			loadSeq(bloomMap, m_kmerSize, m_hashNum, rec.seq, value);
			value++;
			count++;
		}
	}
	if (verbose) {
		std::cerr << "Loaded " << count << " reads from `"
			<< filename << "` into bloom filter\n";
	}

	 bloomMap.storeFilter(filename + ".bm");
}



BloomMapGenerator::~BloomMapGenerator() {
}

