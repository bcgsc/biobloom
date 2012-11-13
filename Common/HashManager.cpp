/*
 * HashManager.cpp
 *
 *  Created on: Aug 10, 2012
 *      Author: cjustin
 */

#include "HashManager.h"

/* Adds a hash function reference to a vector and a seed value to a vector
 * in the same order so they can be executed sequentially
 */
void HashManager::addHashFunction(string const &functionName,
		size_t seedValue) {
	//todo: use ENUM instead of string matching code style
	if (functionName == "CityHash64") {
		hashFuncs.push_back(new CityHash64(seedValue));
	} else {
		cerr << "HashManager: " << functionName
				<< " hash function name not recognized.";
		exit(0);
	}
	tempHashValues.resize(hashFuncs.size());
}

const vector<size_t> &HashManager::multiHash(string const &kmer){
	for (unsigned short i = 0; i < hashFuncs.size(); ++i) {
		//dereference functor call
		tempHashValues[i] = (*hashFuncs[i])(kmer);
	}
	return tempHashValues;
}

HashManager::~HashManager() {
}

