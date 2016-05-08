/*
 * SpacedSeedIndex.h
 *
 *  Created on: Apr 13, 2015
 *      Author: cjustin
 */

#ifndef SPACEDSEEDINDEX_H_
#define SPACEDSEEDINDEX_H_

#include <sstream>
#include <iostream>
#include <fstream>
#include <omp.h>
#include <math.h>
#include <limits>
#include "boost/unordered/unordered_map.hpp"
#include "boost/shared_ptr.hpp"
#include <vector>
#include <stdint.h>
#include "Options.h"

using namespace std;

static const uint8_t b2f[256] = {
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, //0
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, //1
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, //2
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, //3
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
    0xFF, 0x00, 0xFF, 0x01, 0xFF, 0xFF, 0xFF, 0x02, //4   'A' 'C' 'G'
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
    0xFF, 0xFF, 0xFF, 0xFF, 0x03, 0xFF, 0xFF, 0xFF, //5   'T'
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
    0xFF, 0x00, 0xFF, 0x01, 0xFF, 0xFF, 0xFF, 0x02, //6   'a' 'c' 'g'
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
    0xFF, 0xFF, 0xFF, 0xFF, 0x03, 0xFF, 0xFF, 0xFF, //7   't'
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, //8
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, //9
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, //10
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, //11
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, //12
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, //13
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, //14
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, //15
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF
};

static const uint8_t b2r[256] = {
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, //0
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, //1
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, //2
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, //3
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
    0xFF, 0x03, 0xFF, 0x02, 0xFF, 0xFF, 0xFF, 0x01, //4   'A' 'C' 'G'
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
    0xFF, 0xFF, 0xFF, 0xFF, 0x00, 0xFF, 0xFF, 0xFF, //5   'T'
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
    0xFF, 0x03, 0xFF, 0x02, 0xFF, 0xFF, 0xFF, 0x01, //6   'a' 'c' 'g'
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
    0xFF, 0xFF, 0xFF, 0xFF, 0x00, 0xFF, 0xFF, 0xFF, //7   't'
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, //8
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, //9
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, //10
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, //11
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, //12
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, //13		cerr << m_indexSize << endl;
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, //14
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, //		cerr << m_indexSize << endl;15
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF
};

struct entry {
	ID *id;
	uint16_t size;
	uint16_t count;
};

//static const unsigned EMPTY_ID = numeric_limits<unsigned>::max();
static const unsigned DELETED_ENTRY = numeric_limits<uint16_t>::max();
typedef uint64_t SSeed;

template<typename ID>
class SpacedSeedIndex {
public:
	SpacedSeedIndex(const string &seed, unsigned numReads) :
			m_totalSSCount(0), m_removedCount(0), m_readLengths(
					vector<size_t>(numReads,0)), m_seed(seed) {
		unsigned ssWeight = 0;
		for (unsigned i = 0; i < seed.length(); i++)
			if (seed[i] == '1')
				ssWeight++;

		m_indexSize = (SSeed) 1 << 2 * ssWeight;
		m_indTable = (entry *) malloc(m_indexSize * sizeof(entry));
		initTable(m_indTable);
	}

	inline void insert(ID id, const string& seq, uint64_t *occ){
        m_readLengths[id] += seq.length();
        loadSeq(id, seq, occ, m_indTable);
	}

	const inline entry &at(unsigned i) const
	{
		return m_indTable[i];
	}

	/*
	 * Computes mean and stdev of counts
	 */
	//TODO: Check for precision errors (use log()?)
	inline void calCountStats(double &mean, double &stdev) const{
		mean = double(m_totalSSCount*2)/double(m_indexSize);
		double squaredSum = 0;
		for (size_t i = 0; i < m_indexSize; i++) {
			double val = double(this->at(i).count)-mean;
			squaredSum += val * val;
		}
		stdev = sqrt(squaredSum/double(m_indexSize -1));
	}

	inline size_t getTotalCounts() const {
		return m_totalSSCount;
	}

	/*
	 * topKey is the key that is key with the highest multiplicity
	 */
	inline unsigned removeSpacedSeeds(unsigned countThresh, SSeed &topKey){
		unsigned removalCount = 0;
		unsigned topCount = 0;
		for (SSeed i = 0; i < m_indexSize; i++) {
			if (m_indTable[i].count > countThresh
					&& m_indTable[i].count != DELETED_ENTRY) {
				if(topCount < m_indTable[i].count){
					topCount = m_indTable[i].count;
					topKey = i;
				}
				m_indTable[i].count = DELETED_ENTRY;
				m_indTable[i].size = 0;
				free(m_indTable[i].id);
				++removalCount;
			}
		}
		return removalCount;
	}

	/*
	 * Total number of spaced seeds (including dupes) added to index
	 */
	inline unsigned getCounts() const{
		return m_totalSSCount;
	}

	/*
	 * Return size of index
	 */
	inline unsigned size() const{
		return m_indexSize;
	}

	inline size_t getReadLength(ID readID) const{
		return m_readLengths[readID];
	}

	virtual ~SpacedSeedIndex(){
	    // destroy locks
	    for(unsigned ii = 0; ii < m_lockSize; ii++) omp_destroy_lock(&m_locks[ii]);
	    	free(m_locks);
		//free internal arrays
		for (size_t i = 0; i < m_indexSize; i++) {
			if(m_indTable[i].count != 0 && m_indTable[i].count != DELETED_ENTRY)
			{
				free(m_indTable[i].id);
			}
		}
		free(m_indTable);
	}
private:
	SSeed m_indexSize;

	entry *m_indTable;
	typedef boost::shared_ptr< vector<SSeed> > SeedSet;
	size_t m_totalSSCount;

	//TODO: MAKE SURE VALUE IS ADDED UP CORRECTLY (CURRENTLY DOES NOT COUNT REDUNDANT ENTRIES!)
	size_t m_removedCount;

	vector<size_t> m_readLengths;

	const string &m_seed;

	unsigned m_lockSize;
	omp_lock_t *m_locks;

	//returns internal index location
	inline SSeed hashInsert(SSeed j, ID element, entry *T){
	    if (T[j].id == NULL) {
	        T[j].size = 1;
	        T[j].id = (ID *) malloc(T[j].size * sizeof(ID));
	        T[j].id[0] = element;
	        return 0;
	    }
	    else {
	        ++T[j].size;
	        T[j].id = (ID *) realloc(T[j].id, T[j].size * sizeof(ID));
	        T[j].id[T[j].size-1] = element;
	        return T[j].size-1;
	    }
	}

//	//returns internal index location
//	inline SSeed rhashInsert(SSeed j, ID element, entry *T){
//	    if (T[j].rid == NULL) {
//	        T[j].rsize = 1;
//	        T[j].rid = (ID *) malloc(T[j].rsize * sizeof(ID));
//	        T[j].rid[0] = element;
//	        return 0;
//	    }
//	    else {
//	        ++T[j].rsize;
//	        T[j].rid = (ID *) realloc(T[j].rid, T[j].rsize * sizeof(ID));
//	        T[j].rid[T[j].rsize-1] = element;
//	        return T[j].rsize-1;
//	    }
//	}

	inline void loadSeq(const ID readId, const string &faSeq, uint64_t *occ, entry *indTable) {
		if (faSeq.length() < m_seed.length())
	        return;
		unsigned length = faSeq.length()-m_seed.length()+1;

		for (unsigned i=0; i < length; i++) {
			SSeed fSeed = 0, rSeed = 0;
			bool nonACGT = false;
			for (unsigned j = 0; j < m_seed.length(); j++) {
				if (m_seed[j] == '1') {
					uint8_t fChar = b2f[(uint8_t) faSeq[j + i]];
					uint8_t rChar = b2r[(uint8_t) faSeq[m_seed.length() - 1
							- j + i]];
					if (fChar == 0xFF || rChar == 0xFF) {
						nonACGT = true;
						break;
					}
					fSeed = fSeed << 2 | fChar;
					rSeed = rSeed << 2 | rChar;
				}
			}
	        if(nonACGT) continue;

	        SSeed seed = rSeed < fSeed ? rSeed : fSeed ;

			if ((occ[seed / 64] & ((uint64_t) 1 << (63 - seed % 64))) == 0) {
				omp_set_lock(&m_locks[(uint16_t) seed]);
//				if (indTable[seed].count != DELETED_ENTRY) {
//					if (indTable[seed].count > opt::maxCount) {
//						indTable[seed].count = DELETED_ENTRY;
//						indTable[seed].size = 0;
//						free(indTable[seed].id);
//#pragma omp atomic
//						m_removedCount += opt::maxCount;
//					} else {
						hashInsert(seed, readId, indTable);
						++indTable[seed].count;
//					}
//				}
				omp_unset_lock(&m_locks[(uint16_t) seed]);
#pragma omp atomic
				++m_totalSSCount;
				occ[seed / 64] |= ((uint64_t) 1 << (63 - seed % 64));
			}
		}
	}

	inline void initTable(entry *indTable){
	    for (unsigned i=0; i < m_indexSize; i++) { // replacing with calloc??
			indTable[i].id = NULL;
	        indTable[i].size = 0;
	        indTable[i].count = 0;
		}
	    //TODO do no hardcode number of lock
	    unsigned m_lockSize = 65536;
	    // parallelize hash loading using omp_lock_t
//	    unsigned lockSize = opt::lockSize < m_indexSize ? opt::lockSize : m_indexSize;
	    m_locks = (omp_lock_t *) malloc(m_lockSize * sizeof(omp_lock_t));
	    for(unsigned ii = 0; ii < m_lockSize; ii++)
	        omp_init_lock(&m_locks[ii]);
	}

	inline void getindTable(const char *aName, entry *indTable) {

	    ifstream rFile(aName);
	    unsigned readId = 0;
	    #pragma omp parallel
	    {
	    	//table for flagging occupancy of each readID/thread being inserted for each spaced seed
	    	//faster than using table directly
	        uint64_t *occ = (uint64_t *) malloc(((m_indexSize + 63)/64) * sizeof(uint64_t));
	        unsigned lreadId=0;
	        for (;;) {
	            bool good;
	            string faHead, faSeq;
	            #pragma omp critical(rFile)
	            {
	                good=getline(rFile, faHead);
	                good=getline(rFile, faSeq);
	                lreadId = readId;
	                readId++;
	            }
	            if(!good) break;
	            //check for reads shorter than k-merSize
				if (faSeq.length() < m_seed.size()) {
					cerr << "ignoring short read; Length: " << faSeq.length() << " ID: " << faHead << endl;
					continue;
				}
	            for (unsigned i=0; i < (m_indexSize + 63)/64; ++i) occ[i]=0;
	            m_readLengths[lreadId] = faSeq.length();
	            loadSeq(lreadId, faSeq, occ, indTable);
#pragma omp atomic
	            m_totalSSCount += faSeq.length() - m_seed.size() + 1;
	        }
	        free(occ);
	    }
	    rFile.close();

	    unsigned fullBin=0;
	    for (unsigned i=0; i < m_indexSize; i++)
	        if (indTable[i].size!=0) ++fullBin;
	    cerr << "Entries in index table:" << fullBin << "\n";
	}
};

#endif /* SPACEDSEEDINDEX_H_ */
