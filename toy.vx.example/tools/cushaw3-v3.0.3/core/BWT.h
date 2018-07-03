#ifndef _BWT_H
#define _BWT_H
#include "Macros.h"

class BWT
{
public:
	BWT(const char* bwtFile);
	~BWT();

	//general inline functions
	inline uint32_t getBwtSeqLength() {
		return _bwtSeqLength;
	}
	inline uint32_t getDollarPos() {
		return _bwtDollar;
	}
	inline uint32_t* getCCounts() {
		return _bwtCCounts;
	}
	inline uint32_t* getBWTPtr() {
		return _bwtPtr;
	}
	inline uint32_t getBWTSize() {
		return _bwtSize;
	}

	//calculate the cumulative occurrence of a base
	uint32_t bwtOcc(uint8_t base, uint32_t pos);

	//calculate the marked position and its offset to the real position
	inline uint32_t bwtGetMarkedPos(uint32_t saFactor, uint32_t pos,
			uint32_t& mapOff) {
		uint32_t off = 0;
		while (pos % saFactor != 0) {
			++off;
			pos = bwtLFCM(pos);
		}
		mapOff = off;
		return pos;
	}

private:
	/*private member variables*/
	uint32_t _bwtSeqLength; //the number of bases of the compressed sequence
	uint32_t _bwtDollar; //the position of $ symbol in the BWT, which has been excluded
	uint32_t _bwtCCounts[BWT_NUM_OCC]; //Cumulative counts for each base
	uint32_t _bwtSize; //size of the BWT-based sequence
	uint32_t* _bwtPtr; //burrows-wheeler transform for the sequence

private:
	/*private member functions*/
	//get the character at the specified position of the BWT
	inline uint8_t bwtGetBase(uint32_t bwtPos) {
		uint32_t shift;
		uint32_t* baseAddr;
		baseAddr = _bwtPtr
				+ (bwtPos >> BWT_OCC_INTERVAL_SHIFT) * BWT_OCC_PTR_OFFSET;
		baseAddr += BWT_NUM_NUCLEOTIDE;
		baseAddr += ((bwtPos & BWT_OCC_INTERVAL_MASK) >> 4);
		shift = (~bwtPos & 0x0F) << 1;

		return (*baseAddr >> shift) & 3;
	}
	//perform one iteration of "last-to-first column mapping"
	inline uint32_t bwtLFCM(uint32_t bwtPos) {
		uint8_t base;
		if (bwtPos == _bwtDollar)
			return 0;
		if (bwtPos < _bwtDollar) {
			base = bwtGetBase(bwtPos);
			return _bwtCCounts[base] + bwtOcc(base, bwtPos);
		}
		base = bwtGetBase(bwtPos - 1);
		return _bwtCCounts[base] + bwtOcc(base, bwtPos);
	}

	//get the nearest cumulative occurrence of a base
	inline uint32_t* bwtInterleavedOcc(uint32_t pos, uint8_t base) {
		return _bwtPtr + (pos >> BWT_OCC_INTERVAL_SHIFT) * BWT_OCC_PTR_OFFSET;
	}

	/*We use population count to calculate the number of occurrences of a character*/
	inline int bwtOccAux(uint64_t y, int c) {
		// reduce nucleotide counting to bits counting
		y = ((c & 2) ? y : ~y) >> 1 & ((c & 1) ? y : ~y)
				& 0x5555555555555555ull;
		// count the number of 1s in y
		y = (y & 0x3333333333333333ull) + (y >> 2 & 0x3333333333333333ull);

		return ((y + (y >> 4)) & 0xf0f0f0f0f0f0f0full) * 0x101010101010101ull
				>> 56;
	}

	//load the BWT from the file
	void loadFromFile(const char* fileName);

	/*friend classes*/
	friend class MemEngine;
};
#endif
