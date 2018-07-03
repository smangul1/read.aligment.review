/*
 * CigarAlign.h
 *
 *  Created on: Dec 30, 2011
 *      Author: yongchao
 */

#ifndef CIGARALIGN_H_
#define CIGARALIGN_H_
#include "Macros.h"
#include "Utils.h"
#include "Sequence.h"

/*alignment operations*/
#define ALIGN_OP_M	0		/*means substitution*/
#define ALIGN_OP_I	1		/*means insertion at the first sequence*/
#define ALIGN_OP_D	2		/*means deletion at the first sequence*/
#define ALIGN_OP_S	3		/*means the part out of the local alignment*/

class CigarAlign
{
public:
	CigarAlign(Sequence& seq, int match);
	CigarAlign(int32_t alignType, int32_t alignLength, int32_t start1, int32_t end1,
			int32_t start2, int32_t end2, int32_t numMismatches,
			int32_t numIndels, int32_t numGaps, uint32_t* cigar, int32_t ncigar,
			int32_t alignScore, bool rev);
	~CigarAlign();

	inline uint32_t* getCigar() {
		return _cigar;
	}
	inline int32_t getCigarLength() {
		return _ncigar;
	}
	inline int32_t getAlignType()
	{
		return _alignType;
	}
	inline int32_t getEditDistance() {
		return _numMismatches + _numIndels;
	}
	inline int32_t getNumMismatches() {
		return _numMismatches;
	}
	inline int32_t getNumIndels() {
		return _numIndels;
	}
	inline int32_t getNumGaps() {
		return _numGaps;
	}
	inline float getIdentity() {
		/*return value in the range [0, 100)*/
		return ((float) (_alignLength - _numMismatches - _numIndels) * 100)
				/ _alignLength;
	}
	inline int32_t getAlignLength() {
		return _alignLength;
	}
	inline int32_t getAlignScore() {
		return _alignScore;
	}
	/*get the alignment information*/
	inline int32_t getNumBases1() {
		return _end1 - _start1 + 1;
	}
	inline int32_t getNumBases2() {
		return _end2 - _start2 + 1;
	}
	inline int32_t getStart() {
		return _start1;
	}
	inline int32_t getEnd() {
		return _end1;
	}
	inline int32_t getMateStart() {
		return _start2;
	}
	inline int32_t getMateEnd() {
		return _end2;
	}
	inline void getRegion(int32_t& start1, int32_t& start2, int32_t& end1,
			int32_t& end2) {
		start1 = _start1;
		start2 = _start2;
		end1 = _end1;
		end2 = _end2;
	}
	inline void cigarOut(FILE* file) {
		uint32_t cigar;
		if (file == NULL) {
			Utils::exit("Invalid file pointer in function %s line %d\n",
					__FUNCTION__, __LINE__);
		}
		for (int32_t i = 0; i < _ncigar; ++i) {
			cigar = _cigar[i];
			/*print the cigar to the file*/
			fprintf(file, "%d%c", cigar >> 2, _alignOpName[cigar & 3]);
		}
	}
	inline void cigarOut(gzFile file) {
		uint32_t cigar;
		if (file == NULL) {
			Utils::exit("Invalid file pointer in function %s line %d\n",
					__FUNCTION__, __LINE__);
		}
		for (int32_t i = 0; i < _ncigar; ++i) {
			cigar = _cigar[i];
			/*print the cigar to the file*/
			gzprintf(file, "%d%c", cigar >> 2, _alignOpName[cigar & 3]);
		}
	}
	/*extend the cigar to the full length of the sequence*/
	int extendCigar(int32_t length);

private:
	int32_t _start1, _end1; /*its self*/
	int32_t _start2, _end2; /*its mate*/
	int32_t _alignLength; /*alignment length*/
	int32_t _alignScore;
	int32_t _numMismatches; /*number of mismatches*/
	int32_t _numIndels; /*number of indels*/
	int32_t _numGaps; /*number of gaps*/

	uint32_t* _cigar;
	int32_t _ncigar;
	int32_t _alignType;

	/*reverse operations for cigar*/
	static const uint8_t _alignRcOp[4];
	static const uint8_t _alignOpName[4];
};

#endif /* CigarAlign_H_ */
