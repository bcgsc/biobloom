/*
 * BloomMapGenerator.cpp
 *
 *  Created on: Mar 17, 2016
 *      Author: cjustin
 */

#include "BloomMapGenerator.h"
#include "DataLayer/FastaReader.h"

/*
 * If size is set to 0 (or not set) a filter size will be estimated based
 * on the number of k-mers in the file
 */
BloomMapGenerator::BloomMapGenerator(vector<string> const &filenames,
		unsigned kmerSize, unsigned hashNum, size_t numElements = 0) :
		m_kmerSize(kmerSize), m_hashNum(hashNum),  m_expectedEntries(
				numElements), m_totalEntries(0), m_fileNames(filenames) {
	//Instantiate dense hash map
	m_headerIDs.set_empty_key(opt::EMPTY);

	//estimate number of k-mers
	if (numElements == 0) {
		//assume each file one line per file
#pragma omp parallel
		for (vector<string>::const_iterator it = m_fileNames.begin();
				it != m_fileNames.end(); ++it) {
			cerr << "Counting total k-mers in " << *it << endl;
			FastaReader sequence(it->c_str(), FastaReader::NO_FOLD_CASE);
			for (FastqRecord rec;;) {
				bool good;
#pragma omp critical(sequence)
				{
					good = sequence >> rec;
				}
				if (good) {
					m_expectedEntries += rec.seq.length() - m_kmerSize;
				} else
					break;
			}
		}
	}
	cerr << "Expected number of elements: " << m_expectedEntries << endl;
}

/* Generate the bloom filter to the output filename
 */
void BloomMapGenerator::generate(const string &filename, double fpr) {
	//init bloom map
	BloomMap<ID> bloomMap(m_expectedEntries, fpr, m_hashNum, m_kmerSize);
	ID value = 0;

	//for each file
	for (vector<string>::const_iterator it = m_fileNames.begin();
			it != m_fileNames.end(); ++it) {
//		if (opt::verbose)
			std::cerr << "Reading `" << *it << "'...\n";
		//read file
		FastaReader sequence(it->c_str(), FastaReader::NO_FOLD_CASE);
		for (FastqRecord rec;;) {
			bool good;
			{
				good = sequence >> rec;
				value++;
			}
			if (good) {
				//k-merize with rolling hash insert into bloom map
				loadSeq(bloomMap, rec.seq, value);
				//assign header to unique ID
				m_headerIDs[value] = rec.id;
			} else
				break;
		}
	}

	//save filter
	bloomMap.storeFilter(filename);
	writeIDs("out.txt", m_headerIDs);
}



BloomMapGenerator::~BloomMapGenerator() {
}

