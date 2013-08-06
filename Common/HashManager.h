/*
 * HashManager.h
 *
 *  Created on: Aug 10, 2012
 *      Author: cjustin
 */
//@TODO: experiment with hash concepts by Adam Kirsch and Michael Mitzenmacher in Building a Better Bloom Filter

#ifndef HASHMANAGER_H_
#define HASHMANAGER_H_
#include <vector>
#include "city.h"

using namespace std;

static inline vector<size_t> multiHash(const char* kmer,size_t kmerSize, size_t num) {
	vector<size_t> tempHashValues(num);
	for (size_t i = 0; i < num; ++i) {
		tempHashValues[i] = CityHash64WithSeed(kmer, kmerSize, i );
	}
	return tempHashValues;
}

#endif /* HASHMANAGER_H_ */
