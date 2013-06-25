/*
 * HashManager.h
 *
 *  Created on: Aug 10, 2012
 *      Author: cjustin
 */

#ifndef HASHMANAGER_H_
#define HASHMANAGER_H_
#include <vector>
#include <string>
#include <iostream>
#include "boost/dynamic_bitset.hpp"
#include "city.h"

using namespace std;

class HashManager {
public:
	void addHashFunction(string const &functionName, size_t seedValue);
	vector<size_t> &multiHash(string const &kmer) const;
	//for possibly binary version of DNA representation
	//vector<size_t> multiHash(boost::dynamic_bitset kmer);

	virtual ~HashManager();
private:
	//the base hash class
	class Hash {
	public:
		Hash() {
		}
		Hash(size_t seedValue) {
			seedVal = seedValue;
		}
		virtual const size_t operator()(string const &kmer) const=0;
		virtual ~Hash() {
		}
	protected:
		size_t seedVal;
	};

	//functor for storing hash functions to be run
	class CityHash64: public Hash {
	public:
		CityHash64(size_t seedValue) {
			seedVal = seedValue;
		}
		const size_t operator()(string const &kmer) const {
			return CityHash64WithSeed(kmer.c_str(), kmer.length(), seedVal);
		}
		virtual ~CityHash64() {
		}
	};
	vector<Hash*> hashFuncs;
};

#endif /* HASHMANAGER_H_ */
