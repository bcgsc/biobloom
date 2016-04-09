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
		unsigned kmerSize, size_t numElements = 0) :
		m_kmerSize(kmerSize), m_expectedEntries(numElements), m_totalEntries(0), m_fileNames(
				filenames), m_collisionCount(0){
	//Instantiate dense hash map
	m_headerIDs.set_empty_key(opt::EMPTY);
//	m_collisionIn.set_empty_key(opt::EMPTY);
//	m_collisionOut.set_empty_key(opt::EMPTY << ID_BITS | opt::EMPTY);

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
//					++m_currentID;
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
void BloomMapGenerator::generate(const string &filePrefix, double fpr) {
	//init bloom map
	BloomMapSS<ID> bloomMap(m_expectedEntries, fpr, opt::sseeds);
	ID value = 0;
	size_t uniqueCount = 0;

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
//				cerr << value << endl;;
				//k-merize with rolling hash insert into bloom map
				uniqueCount += loadSeq(bloomMap, rec.seq, value);
				//assign header to unique ID
				m_headerIDs[value] = rec.id;
//				PairID pairID = value << 16 | value;
//				m_collisionIn[value] = pairID;
//				m_collisionOut[pairID] = value;
			} else
				break;
		}
	}
	cerr << "Collision Rate: " << double(m_collisionCount)/double(m_expectedEntries*opt::sseeds.size()) << endl;
	bloomMap.setUnique(uniqueCount);
	cerr << "Compressing Bloom Map" << endl;
	BloomMapSSBitVec<ID> bloomMapBV(bloomMap);

	//save filter
	bloomMap.storeFilter(filePrefix + "_uncompressed.bf");
	bloomMapBV.store(filePrefix + ".bf");
	writeIDs(filePrefix + "_ids.txt", m_headerIDs);
}



BloomMapGenerator::~BloomMapGenerator() {
}

