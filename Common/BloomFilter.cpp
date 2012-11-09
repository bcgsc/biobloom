/*
 * BloomFilter.cpp
 *
 *  Created on: Aug 10, 2012
 *      Author: cjustin
 */

#include "BloomFilter.h"
#include <vector>
#include <fstream>
#include <iostream>

//todo: Performance optimization potential: make size of bitset a power of 2
/* De novo filter constructor.
 * precondition: filterSize must be a multiple of 64
 */
BloomFilter::BloomFilter(size_t const &filterSize, HashManager const &hashFns) :
		multiHasher(hashFns) {
	if (filterSize % 64 != 0) {
		cerr << "ERROR: Filter Size \"" << filterSize
				<< "\" is not a multiple of 64." << endl;
		exit(1);
	}
	filter = boost::dynamic_bitset<>(filterSize);
}

/*
 * Loads the filter (file is a .bf file) from path specified
 */
//todo: OptPot: add decompression/compression support for loading speed
//todo: OptPot: change to forked process for loading speed
BloomFilter::BloomFilter(string const &filterFilePath,
		HashManager const &hashFns) :
		multiHasher(hashFns) {
	vector<boost::dynamic_bitset<>::block_type> filterBlocks;

	//load in blocks to a vector
	ifstream binaryFile(filterFilePath.c_str(), ios::binary);
	if (binaryFile.is_open()) {
		while (binaryFile.good()) {
			char buffer[sizeof(boost::dynamic_bitset<>::block_type)];
			binaryFile.read(buffer,
					sizeof(boost::dynamic_bitset<>::block_type));
			boost::dynamic_bitset<>::block_type block;
			memcpy(&block, buffer, sizeof(block));
			filterBlocks.push_back(block);
		}
		//remove last block which is eof;
		filterBlocks.pop_back();
	} else {
		cerr << "file \"" << filterFilePath << "\" could not be read." << endl;
		binaryFile.close();
		exit(1);
	}
	binaryFile.close();

	//load vector into filter
	//constructor builds bitset size based on number of blocks (multiple of 64)
	filter = boost::dynamic_bitset<>(filterBlocks.begin(), filterBlocks.end());
}

void BloomFilter::insert(string const &kmer) {
	vector<size_t> values = multiHasher.multiHash(kmer);
	//iterates through hashed values adding it to the filter
	for (vector<size_t>::iterator it = values.begin(); it != values.end();
			++it) {
		filter.set(*it % filter.size());
	}
}

/*
 * Stores the filter as a binary file to the path specified
 */
void BloomFilter::storeFilter(string const &filterFilePath) const {
	vector<boost::dynamic_bitset<>::block_type> filterBlocks(
			filter.num_blocks());

	//populate vector blocks
	unsigned long result;
	boost::to_block_range(filter, filterBlocks.begin());

	ofstream myFile(filterFilePath.c_str(), ios::out | ios::binary);

	//todo: Optimization potential: concatenate blocks (after conversion to char*) before writing to increase writing speed
	//todo: OptPot: send data to a compression stream to decrease bytes needed to be written

	//write out each block
	for (vector<boost::dynamic_bitset<>::block_type>::iterator it =
			filterBlocks.begin(); it != filterBlocks.end(); ++it) {

		//retrieves block and converts it to a char*
		myFile.write(reinterpret_cast<char*>(&*it),
				sizeof(boost::dynamic_bitset<>::block_type));
	}
	myFile.close();
}

const bool BloomFilter::contains(string const &kmer) {
	return contains(multiHasher.multiHash(kmer));
}

/*
 * Accepts a list of precomputed hash values. Faster than rehashing each time.
 */
const bool BloomFilter::contains(vector<size_t> const &values) {
	for (vector<size_t>::const_iterator it = values.begin(); it != values.end();
			++it) {
		if (!filter.test(*it % filter.size())) {
			return false;
		}
	}
	return true;
}

BloomFilter::~BloomFilter() {
// TODO Auto-generated destructor stub
}

