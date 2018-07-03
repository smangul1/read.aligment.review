/*
 * CigarAlign.cpp
 *
 *  Created on: Dec 30, 2011
 *      Author: yongchao
 */
#include "CigarAlign.h"
#include "Utils.h"

const uint8_t CigarAlign::_alignRcOp[4] = { ALIGN_OP_M, ALIGN_OP_D, ALIGN_OP_I,
		ALIGN_OP_S };
const uint8_t CigarAlign::_alignOpName[4] = { 'M', 'I', 'D', 'S' };

CigarAlign::CigarAlign(Sequence& seq, int match) {
	/*this constructs an exact match for sequence*/
	_start1 = 0;
	_end1 = seq._length - 1;
	_start2 = _start1;
	_end2 = _end1;
	_alignLength = seq._length;
	_alignScore = seq._length * match;
	_alignType = REALIGN_TYPE_LOCAL;
	_numMismatches = _numIndels = _numGaps = 0;
	_ncigar = 1;
	_cigar = new uint32_t[_ncigar];
	_cigar[0] = (seq._length << 2) | ALIGN_OP_M;
}
CigarAlign::CigarAlign(int32_t alignType, int32_t alignLength, int32_t start1, int32_t end1,
		int32_t start2, int32_t end2, int32_t numMismatches, int32_t numIndels,
		int32_t numGaps, uint32_t* cigar, int32_t ncigar, int32_t alignScore,
		bool isfirst) {
	/*the part in the local alignment for the sequence itself*/
	_start1 = start1;
	_end1 = end1;
	/*the part in the local alignment for the other sequence*/
	_start2 = start2;
	_end2 = end2;

	/*the local alignment length*/
	_alignLength = alignLength;
	_alignScore = alignScore;
	_alignType = alignType;

	/*save the number of errors*/
	_numMismatches = numMismatches;
	_numIndels = numIndels;
	_numGaps = numGaps;

	/*save the cigar CigarAlign*/
	_ncigar = ncigar;
	_cigar = new uint32_t[ncigar];
	if (_cigar == NULL) {
		Utils::exit("Memory allocation failed in function %s line %d\n",
				__FUNCTION__, __LINE__);
	}

	/*copy the data*/
	if (isfirst == true) {
		for (int32_t i = 0, j = ncigar - 1; i < ncigar; ++i, --j) {
			_cigar[i] = cigar[j];
		}
	} else {
		for (int32_t i = 0, j = ncigar - 1; i < ncigar; ++i, --j) {
			uint32_t tcigar = cigar[j];
			_cigar[i] = (tcigar & 0xFFFC) | _alignRcOp[tcigar & 3];
		}
	}
}

CigarAlign::~CigarAlign() {
	delete[] _cigar;
}
int CigarAlign::extendCigar(int32_t length) {
	uint32_t* buffer;

	if (_start1 < 0 || _end1 >= length) {
		Utils::log("BUG: some errors occurred in the banded local alignment (%d %d %d)\n", _start1, _end1, length);
		return 0;
	}
	if (_start1 == 0 && _end1 == length - 1) {
		return 0;
	}

	buffer = new uint32_t[_ncigar + 2];
	if (buffer == NULL) {
		Utils::exit("Memory allocation failed in function %s line %d\n",
				__FUNCTION__, __LINE__);
	}
	//extend to the begining of the sequence without gaps
	if (_start1 > 0) {
		buffer[0] = (_start1 << 2) | ALIGN_OP_S;
		memcpy(buffer + 1, _cigar, _ncigar * sizeof(_cigar[0]));
		_ncigar++;
	} else {
		memcpy(buffer, _cigar, _ncigar * sizeof(_cigar[0]));
	}

	//extend to the end of the sequence without gaps
	if (_end1 < length - 1) {
		buffer[_ncigar++] = ((length - _end1 - 1) << 2) | ALIGN_OP_S;
	}

	/*replace the cigar content*/
	delete[] _cigar;
	_cigar = buffer;

	return 1;
}

