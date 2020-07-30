/*
 * SDust.hpp
 *
 *	Cpp interface with sdust. Assumes position queries will be in order.
 *  Created on: Jul. 29, 2020
 *      Author: cjustin
 */

//TODO: There are unnecessary memory allocations. Please optimize me.
#ifndef COMMON_SDUST_HPP_
#define COMMON_SDUST_HPP_

#include "sdust.h"
#include "Options.h"
#include <utility>

class SDust {
public:
	SDust() :
		m_currPos(0), m_results(0), m_resSize(0) {
	}

	SDust(const string &seq): m_currPos(0) {
		m_results = (uint64_t*)sdust(0, (uint8_t*) seq.c_str(), seq.size(), opt::dustT,
				opt::dustWindow, &m_resSize);
	}

	~SDust() {
		free(m_results);
	}

	void loadSeq(const string &seq) {
		free(m_results);
		m_results = (uint64_t*)sdust(0, (uint8_t*) seq.c_str(), seq.size(), opt::dustT,
				opt::dustWindow, &m_resSize);
	}

	//Assumes queries are in order
	bool isLowComp(unsigned pos) {
		while (m_currPos < m_resSize && pos >= (unsigned) m_results[m_currPos]) {
			++m_currPos;
		}
		if (m_currPos < m_resSize
				&& pos >= (unsigned) (m_results[m_currPos] >> 32)
				&& pos < (unsigned) m_results[m_currPos]) {
			return true;
		}
		return false;
	}

private:
	int m_currPos;
	uint64_t *m_results;
	int m_resSize;
};

#endif /* COMMON_SDUST_HPP_ */
