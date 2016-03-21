/*
 * BloomMapGenerator.cpp
 *
 *  Created on: Mar 17, 2016
 *      Author: cjustin
 */

#include <BloomMapGenerator.h>
#include "DataLayer/FastaReader.h"
#include "bloomfilter/BloomMap.hpp"

/*
 * If size is set to 0 (or not set) a filter size will be estimated based
 * on the number of k-mers in the file
 */
BloomMapGenerator::BloomMapGenerator(vector<string> const &filenames,
		unsigned kmerSize, unsigned hashNum, size_t numElements = 0) :
		m_kmerSize(kmerSize), m_hashNum(hashNum), m_expectedEntries(
				numElements), m_fileNames(filenames) {
	if (numElements == 0) {
		//assume each file one line per file
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

/* Generate the bloom filter to the input filename
 * Returns number of redundant entries
 */
size_t BloomMapGenerator::generate(const string &filename) {
	//init bloom map



}

BloomMapGenerator::~BloomMapGenerator() {
}

