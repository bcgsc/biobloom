#ifndef ABYSS_ROLLING_HASH_H
#define ABYSS_ROLLING_HASH_H 1

#include "rolling.h"
#include <string>
#include <vector>
#include <cassert>

class RollingHash
{
private:

	/**
	 * Determine the canonical hash value, given hash values for
	 * forward and reverse-complement of the same k-mer.
	 */
	inline size_t canonicalHash(size_t hash, size_t rcHash) const
	{
		return (rcHash < hash) ? rcHash : hash;
	}

//	/** compute multiple pseudo-independent hash values using
//	 * a single seed hash value */
//	std::vector<size_t> multiHash(size_t seedHash) const
//	{
//		std::vector<size_t> hashes(m_numHashes);
//		for (unsigned i = 0; i < m_numHashes; i++) {
//			hashes.at(i) = rol(varSeed, i) ^ seedHash;
//		}
//		return hashes;
//	}

	/**
	 * Compute multiple pseudo-independent hash values using spaced seed values
	 * Also computes reverse which complement is found hashes as well
	 * kmerItr is the location in the string you are currently at to use for masking
	 * @return vector of hash values
	 **/
	//TODO: Is the construction of an array really needed? (pass by a reference?) what about thread safety?
	void multiHash(uint64_t rHash, uint64_t fHash,
			const std::string::const_iterator &kmerItr, std::vector<size_t> &hashes)
	{
		for (unsigned i = 0; i < m_ssVals.size(); ++i) {
			hashes[i] = maskValues(rHash, fHash, m_k, m_ssVals.at(i), kmerItr);
		}
	}

public:

//	/**
//	 * Default constructor.
//	 */
//	RollingHash() : m_k(0), m_hash1(0), m_rcHash1(0), m_ssVals(
//					vector < vector<unsigned> >()), m_hashes(m_ssVals.size()) {
//	}

	/**
	 * Constructor. Construct RollingHash object when initial k-mer
	 * is unknown.
	 * @param numHashes number of pseudo-independent hash values to compute
	 * for each k-mer
	 * @param k k-mer length
	 */
	RollingHash(unsigned k,	const vector<vector<unsigned> > &ssVals) :
			m_k(k), m_hash1(0), m_rcHash1(0), m_ssVals(ssVals), m_hashes(m_ssVals.size()) {
	}

//	/**
//	 * Constructor. Construct RollingHash object while specifying
//	 * initial k-mer to be hashed.
//	 * @param kmer initial k-mer for initializing hash value(s)
//	 * @param numHashes number of pseudo-independent hash values to compute
//	 * for each k-mer
//	 * @param k k-mer length
//	 */
//	RollingHash(const std::string& kmer, unsigned numHashes, unsigned k)
//		: m_numHashes(numHashes), m_k(k), m_hash1(0), m_rcHash1(0)
//	{
//		reset(kmer);
//	}

	/** initialize hash values from seq */
	void reset(const std::string::const_iterator &kmerItr)
	{
		/* compute first hash function for k-mer */
		m_hash1 = getFhval(kmerItr, m_k);

		/* compute first hash function for reverse complement
		 * of k-mer */
		m_rcHash1 = getRhval(kmerItr, m_k);

		/* compute hash values */
		multiHash(m_hash1, m_rcHash1, kmerItr, m_hashes);
	}

//	/**
//	 * Compute hash values for a neighbour k-mer on the right,
//	 * without updating internal state.
//	 * @param std::string::const_iterator start position of what base is being removed
//	 * @return vector of hash values for next k-mer
//	 */
//	std::vector<size_t> peekRight(const std::string::const_iterator &kmerItr) const
//	{
//		size_t hash1 = m_hash1;
//		size_t rcHash1 = m_rcHash1;
//
//		/* update first hash function */
//		rollHashesRight(hash1, rcHash1, *kmerItr, *(kmerItr+m_k), m_k);
//
////		/* get seed value for computing rest of the hash functions */
////		size_t seed = canonicalHash(hash1, rcHash1);
//
//		/* compute hash values */
//		return multiHash(hash1, rcHash1, kmerItr);
//	}

//	/**
//	 * Compute hash values for a neighbour k-mer on the left,
//	 * without updating internal state.
//	 * @param charIn leftmost base of k-mer to the left
//	 * @param charOut rightmost base of current k-mer
//	 * @return vector of hash values for next k-mer
//	 */
//	std::vector<size_t> peekLeft(const std::string::const_iterator &kmerItr) const
//	{
//		size_t hash1 = m_hash1;
//		size_t rcHash1 = m_rcHash1;
//
//		/* update first hash function */
//		rollHashesLeft(hash1, rcHash1, *kmerItr, *(kmerItr+m_k), m_k);
//
////		/* get seed value for computing rest of the hash functions */
////		size_t seed = canonicalHash(hash1, rcHash1);
//
//		/* compute hash values */
//		return multiHash(m_hash1, m_rcHash1, kmerItr);
//	}

	/**
	 * Compute hash values for next k-mer to the right and
	 * update internal state.
	 * @param charOut leftmost base of current k-mer
	 * @param charIn rightmost base of next k-mer
	 * @return vector of hash values for next k-mer
	 */
	void rollRight(const std::string::const_iterator &kmerItr)
	{
		/* update first hash function */
		rollHashesRight(m_hash1, m_rcHash1, *kmerItr, *(kmerItr + m_k), m_k);
//		/* get seed value for computing rest of the hash functions */
//		size_t seed = canonicalHash(m_hash1, m_rcHash1);

		/* compute hash values */
		//need to shift right the current position on the string by 1
		multiHash(m_hash1, m_rcHash1, kmerItr + 1, m_hashes);
	}

	/**
	 * Compute hash values for next k-mer to the left and
	 * update internal state.
	 * @param charOut leftmost base of current k-mer
	 * @param charIn rightmost base of next k-mer
	 * @return vector of hash values for next k-mer
	 */
	void rollLeft(const std::string::const_iterator &kmerItr)
	{
		/* update first hash function */
		rollHashesLeft(m_hash1, m_rcHash1, *kmerItr, *(kmerItr + m_k), m_k);

//		/* get seed value for computing rest of the hash functions */
//		size_t seed = canonicalHash(m_hash1, m_rcHash1);

		/* compute hash values */
		//need to shift left the current position on the string by 1
		multiHash(m_hash1, m_rcHash1, kmerItr - 1, m_hashes);
	}

	/** Get hash values for current k-mer */
	const std::vector<size_t>& getHash() const
	{
		assert(!m_hashes.empty());
		return m_hashes;
	}

	/** Equality operator */
	bool operator==(const RollingHash& o) const
	{
		return
			m_ssVals == o.m_ssVals &&
			m_k == o.m_k &&
			m_hashes == o.m_hashes &&
			m_hash1 == o.m_hash1 &&
			m_rcHash1 == o.m_rcHash1;
	}

private:

//	/** number of hash functions to compute at each position */
//	unsigned m_numHashes;
	/** k-mer length */
	unsigned m_k;
	/** value of first hash function for current k-mer */
	size_t m_hash1;
	/** value of first hash function for current k-mer, after
	 * reverse-complementing */
	size_t m_rcHash1;
	//Spaced seed values to use
	const vector< vector<unsigned> > &m_ssVals;
	/** current values for each hash function */
	std::vector<size_t> m_hashes;
};

#endif
