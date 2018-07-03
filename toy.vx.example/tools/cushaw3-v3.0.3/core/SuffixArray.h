/*
 * SuffixArray.h
 *
 *  Created on: Dec 24, 2011
 *      Author: yongchao
 */

#ifndef SUFFIXARRAY_H_
#define SUFFIXARRAY_H_

#include "Macros.h"
#include "BWT.h"

class SuffixArray
{
public:
	SuffixArray(const char* saFileName);
	~SuffixArray();

	/*get the length of the sequence indexed using BWT*/
	inline uint32_t getBwtSeqLength() {
		return _seqLength;
	}

	/*get the factor for reduced suffix array*/
	inline uint32_t getFactor() {
		return _saFactor;
	}
	static uint32_t getFactor(const char* saFileName);

	//get the position from the marked position and its offset
	inline uint32_t getPosition(uint32_t markedPosition, uint32_t markedOff) {
		return markedOff + _saPtr[markedPosition / _saFactor];
	}

	//get the position directly from BWT
	inline uint32_t getPosition(BWT* bwt, uint32_t pos) {
		uint32_t markedPosition;
		uint32_t markedOff;

		//get the marked position and its offset
		markedPosition = bwt->bwtGetMarkedPos(_saFactor, pos, markedOff);

		//get the mapping position
		return markedOff + _saPtr[markedPosition / _saFactor];
	}
private:
	uint32_t _seqLength;
	uint32_t* _saPtr;
	uint32_t _saSize;
	uint32_t _saFactor;
};

#endif /* SUFFIXARRAY_H_ */
