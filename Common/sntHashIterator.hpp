/*
 * HashFunction.hpp
 *
 * 	Created to get around templating issues revolving ntHashIterator and stHashIterator
 *
 *
 *  Created on: Aug. 25, 2020
 *      Author: cjustin
 */

#ifndef COMMON_SNTHASHITERATOR_HPP_
#define COMMON_SNTHASHITERATOR_HPP_
#include "btl_bloomfilter/vendor/ntHashIterator.hpp"
#include <assert.h>

using namespace std;

class sntHashIterator : public ntHashIterator
{
  public:
	/*
	 * Default constructor.
	 */
	sntHashIterator()
	  :
			ntHashIterator() {
	}

	sntHashIterator(const std::string &seq,
			const std::vector<std::vector<unsigned> > &seed, unsigned h,
			unsigned h2, unsigned k, size_t pos = 0) :
			ntHashIterator(seq, h, k, pos) {
		assert(seed.empty());
		assert(h2 == 1);
	}
};

#endif /* COMMON_SNTHASHITERATOR_HPP_ */
