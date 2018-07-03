/*
 * Aligner.cpp
 *
 *  Created on: Dec 29, 2011
 *      Author: yongchao
 We have used and modified some code from the open-source SWIPE algorithm in this file.
 The following details the copyright and licence information of SWIPE.

 SWIPE
 Smith-Waterman database searches with Inter-sequence Parallel Execution

 Copyright (C) 2008-2012 Torbjorn Rognes, University of Oslo, 
 Oslo University Hospital and Sencel Bioinformatics AS

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU Affero General Public License as
 published by the Free Software Foundation, either version 3 of the
 License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Affero General Public License for more details.

 You should have received a copy of the GNU Affero General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.

 Contact: Torbjorn Rognes <torognes@ifi.uio.no>, 
 Department of Informatics, University of Oslo, 
 PO Box 1080 Blindern, NO-0316 Oslo, Norway
 */

#include "Aligner.h"
#include "Utils.h"
#include "SeqFileParser.h"

//#define ALIGN_DEBUG
#define cpuid(f,a,b,c,d) asm("cpuid": "=a" (a), "=b" (b), "=c" (c), "=d" (d) : "a" (f));
Aligner::Aligner(Options* options) {
	/*check the SSSE3 and SSE2*/
	unsigned int a __attribute__ ((unused));
	unsigned int b __attribute__ ((unused));
	unsigned int c, d;

	cpuid(1, a, b, c, d);
	//  printf("cpuid: %08x %08x %08x %08x\n", a, b, c, d);

	/*check sse2*/
	if (!((d >> 26) & 1)) {
		Utils::exit("!!!!Requiring a CPU with SSE2 support!!!!\n");
	}
	/*check SSSE3*/
	if ((c >> 9) & 1) {
		_haveSSSE3 = true;
	} else {
		_haveSSSE3 = false;
		Utils::log(
				"!!!!DO NOT have SSSE3 support on the CPU, resulting in lower speed\n");
	}

	_options = options;
	/*get gap open and extension penalties*/
	_gopen = options->getGapOpen();
	_gext = options->getGapExtend();
	_goe = _gopen + _gext;

	/*score matrix*/
	_match = options->getMatch();
	_mismatch = -options->getMismatch();

	/*for mismatches*/
	for (int i = 0; i < 5; ++i) {
		for (int j = 0; j < 5; ++j) {
			_smat[i][j] = (i == j) ? _match : _mismatch;
		}
	}
	/*for unknown bases*/
	for (int i = 0; i < 5; ++i) {
		_smat[UNKNOWN_BASE][i] = _smat[i][UNKNOWN_BASE] = -1;
	}

	_profile = NULL;
	_profileSize = 0;
	_heBuffer = NULL;
	_heBufferSize = 0;
	_tbtable = NULL;
	_tbtableSize = 0;
	_maxQuerySize = MAX_SEQ_LENGTH;

	/*scoring matrix for SSE2*/
	_score_matrix_7 = new char[32 * 32];
	if (_score_matrix_7 == NULL) {
		Utils::exit("Memory allocation failed in function %s line %d\n",
				__FUNCTION__, __LINE__);
	}
	_score_matrix_16 = new int16_t[32 * 32];
	if (_score_matrix_16 == NULL) {
		Utils::exit("Memory allocation failed in function %s line %d\n",
				__FUNCTION__, __LINE__);
	}
	/*initialize the matrix with _smat*/
	int hi = -100;
	for (int i = 1; i <= 5; i++) {
		for (int j = 1; j <= 5; ++j) {
			char score = _smat[i - 1][j - 1];
			_score_matrix_7[(i << 5) + j] = score;
			_score_matrix_16[(i << 5) + j] = score;
			if (score > hi) {
				hi = score;
			}
		}
	}
	_scorelimit7 = 128 - hi;

	_dprofile = new uint8_t[4 * 16 * 32];
	if (_dprofile == NULL) {
		Utils::exit("Memory allocation failed in function %s line %d\n",
				__FUNCTION__, __LINE__);
	}

	_qtable = new BYTE*[_maxQuerySize];
	if (_qtable == NULL) {
		Utils::exit("Memory allocation failed in function %s line %d\n",
				__FUNCTION__, __LINE__);
	}

	_hearray = new BYTE[_maxQuerySize * 32];
	if (_hearray == NULL) {
		Utils::exit("Memory allocation failed in function %s line %d\n",
				__FUNCTION__, __LINE__);
	}
}
Aligner::Aligner(Options* options, int32_t match, int32_t mismatch, int32_t gopen, int32_t gext)
{
	/*check the SSSE3 and SSE2*/
	unsigned int a __attribute__ ((unused));
	unsigned int b __attribute__ ((unused));
	unsigned int c, d;

	cpuid(1, a, b, c, d);
	//  printf("cpuid: %08x %08x %08x %08x\n", a, b, c, d);

	/*check sse2*/
	if (!((d >> 26) & 1))
	{
		Utils::exit("!!!!Requiring a CPU with SSE2 support!!!!\n");
	}
	/*check SSSE3*/
	if ((c >> 9) & 1)
	{
		_haveSSSE3 = true;
	}
	else
	{
		_haveSSSE3 = false;
		Utils::log("!!!!DO NOT have SSSE3 support on the CPU, resulting in lower speed\n");
	}

	_options = options;
	/*get gap open and extension penalties*/
	_gopen = gopen;
	_gext = gext;
	_goe = _gopen + _gext;

	/*score matrix*/
	_match = match;
	_mismatch = -mismatch;

	/*for mismatches*/
	for (int i = 0; i < 5; ++i)
	{
		for (int j = 0; j < 5; ++j)
		{
			_smat[i][j] = (i == j) ? _match : _mismatch;
		}
	}
	/*for unknown bases*/
	for(int i = 0; i < 5; ++i) {
		_smat[UNKNOWN_BASE][i] = _smat[i][UNKNOWN_BASE] = -1;
	}

	_profile = NULL;
	_profileSize = 0;
	_heBuffer = NULL;
	_heBufferSize = 0;
	_tbtable = NULL;
	_tbtableSize = 0;
	_maxQuerySize = MAX_SEQ_LENGTH;

	/*scoring matrix for SSE2*/
	_score_matrix_7 = new char[32 * 32];
	if (_score_matrix_7 == NULL)
	{
		Utils::exit("Memory allocation failed in function %s line %d\n",
				__FUNCTION__, __LINE__);
	}
	_score_matrix_16 = new int16_t[32 * 32];
	if (_score_matrix_16 == NULL)
	{
		Utils::exit("Memory allocation failed in function %s line %d\n",
				__FUNCTION__, __LINE__);
	}
	/*initialize the matrix with _smat*/
	int hi = -100;
	for (int i = 1; i <= 5; i++)
	{
		for (int j = 1; j <= 5; ++j)
		{
			char score = _smat[i - 1][j - 1];
			_score_matrix_7[(i << 5) + j] = score;
			_score_matrix_16[(i << 5) + j] = score;
			if (score > hi)
			{
				hi = score;
			}
		}
	}
	_scorelimit7 = 128 - hi;

	_dprofile = new uint8_t[4 * 16 * 32];
	if (_dprofile == NULL)
	{
		Utils::exit("Memory allocation failed in function %s line %d\n",
				__FUNCTION__, __LINE__);
	}

	_qtable = new BYTE*[_maxQuerySize];
	if (_qtable == NULL)
	{
		Utils::exit("Memory allocation failed in function %s line %d\n",
				__FUNCTION__, __LINE__);
	}

	_hearray = new BYTE[_maxQuerySize * 32];
	if (_hearray == NULL)
	{
		Utils::exit("Memory allocation failed in function %s line %d\n",
				__FUNCTION__, __LINE__);
	}
}

Aligner::~Aligner() {
	if (_profile) {
		delete[] _profile;
	}
	if (_heBuffer) {
		delete[] _heBuffer;
	}
	if (_tbtable) {
		delete[] _tbtable;
	}
	delete[] _qtable;
	delete[] _dprofile;
	delete[] _score_matrix_16;
	delete[] _score_matrix_7;
	delete[] _hearray;
}
//#define ONLY_8CHANNEL
void Aligner::lalignScore(uint8_t* query, int qlen, uint8_t* sequences,
		int* seqOffsets, int numSeqs, AlignScore* scores) {

	if (qlen > _maxQuerySize) {
		_maxQuerySize = qlen * 2;
		if (_qtable) {
			delete[] _qtable;
		}
		_qtable = new BYTE*[_maxQuerySize];
		if (_qtable == NULL) {
			Utils::exit("Memory allocation failed in function %s line %d\n",
					__FUNCTION__, __LINE__);
		}
		if (_hearray) {
			delete[] _hearray;
		}
		_hearray = new BYTE[_maxQuerySize * 32];
		if (_hearray == NULL) {
			Utils::exit("Memory allocation failed in function %s line %d\n",
					__FUNCTION__, __LINE__);
		}
	}

	/*copy the query*/
	for (int32_t i = 0; i < qlen; ++i) {
		_qtable[i] = _dprofile + 64 * query[i];
	}

	/*perform the search*/
#ifdef ONLY_8CHANNEL
	search16((WORD**) _qtable, _goe, _gext, (WORD*) _score_matrix_16,
			(WORD*) _dprofile, (WORD*) _hearray, qlen, numSeqs, seqOffsets,
			sequences, scores);
#else
	if (search7(_qtable, _goe, _gext, (BYTE*) _score_matrix_7, _dprofile,
			_hearray, qlen, numSeqs, seqOffsets, sequences, scores)
			>= _scorelimit7) {
		search16((WORD**) _qtable, _goe, _gext, (WORD*) _score_matrix_16,
				(WORD*) _dprofile, (WORD*) _hearray, qlen, numSeqs, seqOffsets,
				sequences, scores);
	}
#endif
}
int32_t Aligner::lalignScore(uint8_t* s1, uint8_t* s2, int32_t s1Length,
		int32_t s2Length, int32_t low, int32_t up) {
	/******************************************
	 * refer to the paper: Kun-Mao Chao, William R. Pearson and Webb Miller (1992)
	 * Aligner two sequences within a specified diagonal band".
	 * Comput Appl Biosci, 8(5): 481-487 */

	int32_t band;
	int32_t leftd, rightd;
	int32_t lowrow, hirow;
	int32_t score, h, e, f;
	int32_t itrans, jtrans;
	int8_t* mat;
	int2* heBuffer;

	band = up - low + 1;
	if (band < 1) {
		Utils::log(
				"low > up is unacceptable for banded local CigarAlign (%d > %d)\n",
				low, up);
		return 0;
	}
	/*CigarAlign buffer*/
	heBuffer = new int2[band + 2];
	if (heBuffer == NULL) {
		Utils::exit("Memory allocation failed in function %s line %d\n",
				__FUNCTION__, __LINE__);
	}

	/*initialize the diagonals*/
	if (low > 0) {
		leftd = 1;
	} else if (up < 0) {
		leftd = band;
	} else {
		leftd = 1 - low;
	}
	rightd = band;

	/*initialize the rows*/
	lowrow = max(0, -up); /* start index -1 */
	hirow = min(s1Length, s2Length - low); /* end index */

	/*calculate profile*/
	int8_t* p;
	int32_t psize = (s1Length + 1) * 5;
	if (psize > _profileSize) {
		if (_profile) {
			delete[] _profile;
		}
		_profileSize = psize * 2;
		_profile = new int8_t[_profileSize];
	}
	p = _profile + 5;
	for (int32_t i = 0; i < s1Length; ++i) {
		mat = _smat[s1[i]];
		for (int32_t j = 0; j < 5; ++j) {
			p[j] = mat[j];
		}
		p += 5;
	}
	/*initialize the vectors*/
	for (int32_t j = leftd; j <= rightd; ++j) {
		heBuffer[j].x = 0;
		heBuffer[j].y = -_gopen;
	}

	heBuffer[rightd + 1].x = ALIGN_MIN_SCORE;
	heBuffer[rightd + 1].y = ALIGN_MIN_SCORE;

	heBuffer[leftd - 1].x = ALIGN_MIN_SCORE;
	heBuffer[leftd].y = -_gopen;

	score = 0;
	for (int i = lowrow + 1; i <= hirow; ++i) {
		if (leftd > 1) {
			--leftd;
		}

		if (i > s2Length - up) {
			--rightd;
		}

		/*calculate the E value*/
		mat = _profile + i * 5;
		if ((h = heBuffer[leftd + 1].x - _goe)
				> (e = heBuffer[leftd + 1].y - _gext)) {
			e = h;
		}

		/*convert to the original CigarAlign matrix*/
		if ((itrans = leftd + low + i - 1) > 0) {
			h = heBuffer[leftd].x + mat[s2[itrans - 1]];
		}
		if (e > h) {
			h = e;
		}
		if (h < 0) {
			h = 0;
		}

		f = h - _gopen;
		heBuffer[leftd].x = h;
		heBuffer[leftd].y = e;
		if (h > score) {
			score = h;
		}

		for (int32_t curd = leftd + 1; curd <= rightd; ++curd) {
			if ((h = h - _goe) > (f = f - _gext)) {
				f = h;
			}
			if ((h = heBuffer[curd + 1].x - _goe)
					> (e = heBuffer[curd + 1].y - _gext)) {
				e = h;
			}
			jtrans = curd + low + i - 1;
			h = heBuffer[curd].x + mat[s2[jtrans - 1]];
			if (f > h) {
				h = f;
			}
			if (e > h) {
				h = e;
			}
			if (h < 0) {
				h = 0;
			}
			heBuffer[curd].x = h;
			heBuffer[curd].y = e;
			if (h > score) {
				score = h;
			}
		}
	}
	delete[] heBuffer;

	return score;
}

vector<CigarAlign*> Aligner::lalignPath(uint8_t* s1, uint8_t* s2,
		int32_t s1Length, int32_t s2Length, int32_t low, int32_t up,
		int which) {
	/******************************************
	 * refer to the paper: Kun-Mao Chao, William R. Pearson and Webb Miller (1992)
	 * Aligner two sequences within a specified diagonal band".
	 * Comput Appl Biosci, 8(5): 481-487 */

	int32_t band;
	int32_t leftd, rightd;
	int32_t lowrow, hirow;
	int32_t score, h, e, f;
	int32_t jtrans;
	int8_t* mat;
	uint8_t *table; /*trace-back table*/
	int32_t cibest, cjbest;
	int32_t bandBufferSize;
	int8_t dir = ALIGN_DIR_STOP;
	vector<CigarAlign*> aligns;
	int bufferSize;
	aligns.reserve(2);

	band = up - low + 1;
	if (band < 1) {
		Utils::log(
				"low > up is unacceptable for banded local CigarAlign (%d > %d)\n",
				low, up);
		return aligns;
	}
	bandBufferSize = band + 2;

	/*CigarAlign buffer*/
	if (bandBufferSize > _heBufferSize) {
		if (_heBuffer) {
			delete[] _heBuffer;
		}
		_heBufferSize = bandBufferSize * 2;
		_heBuffer = new int2[_heBufferSize];
		if (_heBuffer == NULL) {
			Utils::exit("Memory allocation failed in function %s line %d\n",
					__FUNCTION__, __LINE__);
		}
	}
	/*trace-back table*/
	bufferSize = bandBufferSize * (s1Length + 1);
	if (bufferSize > _tbtableSize) {
		if (_tbtable) {
			delete[] _tbtable;
		}
		_tbtableSize = bufferSize * 2;
		_tbtable = new uint8_t[_tbtableSize];
		if (_tbtable == NULL) {
			Utils::exit("Memory allocation failed in function %s line %d\n",
					__FUNCTION__, __LINE__);
		}
	}

	/*initialize the diagonals*/
	if (low > 0) {
		leftd = 1;
	} else if (up < 0) {
		leftd = band;
	} else {
		leftd = 1 - low;
	}
	rightd = band;

	/*initialize the rows*/
	lowrow = max(0, -up); /* start index -1 */
	hirow = min(s1Length, s2Length - low); /* end index */

	/*initialize the CigarAlign table*/

	table = _tbtable;
	for (int32_t i = 0; i < bandBufferSize; ++i) {
		table[i] = ALIGN_DIR_STOP;
	}
	table = _tbtable + s1Length * bandBufferSize;
	for (int32_t i = 0; i < bandBufferSize; ++i) {
		table[i] = ALIGN_DIR_STOP;
	}
	table = _tbtable;
	for (int32_t i = 0; i <= s1Length; ++i) {
		*table = ALIGN_DIR_STOP;
		table += bandBufferSize;
	}
	table = _tbtable + bandBufferSize - 1;
	for (int32_t i = 0; i <= s1Length; ++i) {
		*table = ALIGN_DIR_STOP;
		table += bandBufferSize;
	}
	/*compute profile*/
	int8_t* p;
	bufferSize = (s2Length + 1) * 5;
	if (bufferSize > _profileSize) {
		if (_profile) {
			delete[] _profile;
		}
		_profileSize = bufferSize * 2;
		_profile = new int8_t[_profileSize];
		if (_profile == NULL) {
			Utils::exit("Memory allocation failed in function %s line %d\n",
					__FUNCTION__, __LINE__);
		}
	}
	p = _profile;
	for (int32_t i = 0; i < 5; ++i) { /*row*/
		mat = _smat[i];
		p++;
		for (int32_t j = 0; j < s2Length; j++) { /*column*/
			*p++ = mat[s2[j]];
		}
	}

	/*initialize the vectors*/
	for (int32_t i = leftd; i <= rightd; ++i) {
		_heBuffer[i].x = 0;
		_heBuffer[i].y = -_gopen;
	}

	_heBuffer[rightd + 1].x = ALIGN_MIN_SCORE;
	_heBuffer[rightd + 1].y = ALIGN_MIN_SCORE;
	_heBuffer[leftd - 1].x = ALIGN_MIN_SCORE;
	_heBuffer[leftd - 1].y = -_gopen;

	/*start the core CigarAlign loop*/
	score = 0; /*the best score*/
	cibest = cjbest = 0; /*the INVALID coordinate for the best score*/
	for (int i = lowrow + 1; i <= hirow; ++i) {
		if (leftd > 1) {
			--leftd;
		}

		if (i > s2Length - up) {
			--rightd;
		}

		/*calculate the E value*/
		h = _heBuffer[leftd + 1].x - _goe;
		e = _heBuffer[leftd + 1].y - _gext;
		if (h > e) {
			e = h;
		}
		//mat = _smat[s1[i - 1]];
		mat = _profile + s1[i - 1] * (s2Length + 1);
		if ((jtrans = leftd + low + i - 1) > 0) {
			h = _heBuffer[leftd].x + mat[jtrans];
			dir = ALIGN_DIR_DIAGONAL;
		}

		if (e > h) {
			h = e;
			dir = ALIGN_DIR_UP;
		}
		if (h < 0) {
			h = 0;
			dir = ALIGN_DIR_STOP;
		}

		f = h - _gopen;
		_heBuffer[leftd].x = h;
		_heBuffer[leftd].y = e;

		/*save the score and its coordinate*/
		if (h > score) {
			score = h;
			cibest = i;
			cjbest = leftd;
		}

		/*save the trace-back coordinate*/
		table = _tbtable + i * bandBufferSize;
		table[leftd] = dir;

		/*for the other cells in the i-th row*/
		for (int32_t curd = leftd + 1; curd <= rightd; ++curd) {
			if ((h = h - _goe) > (f = f - _gext)) {
				f = h;
			}
			if ((h = _heBuffer[curd + 1].x - _goe)
					> (e = _heBuffer[curd + 1].y - _gext)) {
				e = h;
			}

			dir = ALIGN_DIR_DIAGONAL;
			jtrans = curd + low + i - 1;
			h = _heBuffer[curd].x + mat[jtrans];
			if (f > h) {
				h = f;
				dir = ALIGN_DIR_LEFT;
			}
			if (e > h) {
				h = e;
				dir = ALIGN_DIR_UP;
			}
			if (h < 0) {
				h = 0;
				dir = ALIGN_DIR_STOP;
			}

			_heBuffer[curd].x = h;
			_heBuffer[curd].y = e;
			if (h > score) {
				score = h;
				cibest = i;
				cjbest = curd;
			}
			//Utils::log("%d ", score);
			/*save the trace-back coordinate*/
			table[curd] = dir;
			//Utils::log("i: %d j: %d dir:%d mat: %d h %d score %d\n", i, curd + low + i - 1, dir, mat[s2[jtrans - 1]], h, score);
		}
	}

	//Utils::log("cibest: %d cjbest: %d score: %d\n", cibest, cjbest + low + cibest - 1, score);
	//check the availability of the local CigarAlign*/
	if (cibest == 0 && cjbest == 0) {
		return aligns;
	}

	/*get the alignment*/
	CigarAlign* align;
	switch (which) {
	case 1: /*the first sequence*/
		align = align2cigar(s1, s2, s1Length, s2Length, _tbtable,
						bandBufferSize, cibest, cjbest, score, low, true);
		if(align){
			aligns.push_back(align);
		}
#ifdef ALIGN_DEBUG
		traceback(s1, s2, s1Length, s2Length, _tbtable, bandBufferSize, cibest,
				cjbest, score, low);
#endif
		break;
	case 2: /*the second sequence*/
		align=align2cigar(s1, s2, s1Length, s2Length, _tbtable,
						bandBufferSize, cibest, cjbest, score, low, false);
		if(align){
			aligns.push_back(align);
		}
#ifdef ALIGN_DEBUG
		traceback(s1, s2, s1Length, s2Length, _tbtable, bandBufferSize, cibest,
				cjbest, score, low);
#endif
		break;

	default: /*both the first and the second sequences*/
		align = align2cigar(s1, s2, s1Length, s2Length, _tbtable,
						bandBufferSize, cibest, cjbest, score, low, true);
		if(align){
			aligns.push_back(align);
		}
#ifdef ALIGN_DEBUG
		traceback(s1, s2, s1Length, s2Length, _tbtable, bandBufferSize, cibest,
				cjbest, score, low);
#endif
		align = align2cigar(s1, s2, s1Length, s2Length, _tbtable,
						bandBufferSize, cibest, cjbest, score, low, false);
		if(align){
			aligns.push_back(align);
		}
#ifdef ALIGN_DEBUG
		traceback(s1, s2, s1Length, s2Length, _tbtable, bandBufferSize, cibest,
				cjbest, score, low);
#endif
		break;
	}

	return aligns;
}

/*global alignment backtracking*/
vector<CigarAlign*> Aligner::galignPath(uint8_t* s1, uint8_t* s2,
		int32_t s1Length, int32_t s2Length) {
	int32_t row, col;
	int32_t bufferSize;
	int32_t h, e, f, p, score;
	int32_t cibest, cjbest;
	uint8_t dir;
	bool border;
	uint8_t* pdata;
	int8_t *pmatrix, *profile;
	vector<CigarAlign*> aligns;
	aligns.reserve(2);

	bufferSize = s1Length; /*s1 is always the query and the alignment is done in parallel to the query*/
	if (bufferSize > _heBufferSize) {
		if (_heBuffer) {
			delete[] _heBuffer;
		}
		_heBufferSize = bufferSize * 2; /*double the buffer size*/
		_heBuffer = new int2[_heBufferSize];
		if (_heBuffer == NULL) {
			Utils::exit("Memory allocation failed in function %s line %d\n",
					__FUNCTION__, __LINE__);
		}
	}
	/*trace-back table: each byte stores 4 alignment moves to save memory*/
	_tbRowSize = (s1Length + 3) >> 2; /*number of bytes of each column*/
	bufferSize = _tbRowSize * s2Length; /*number of rows*/
	if (bufferSize > _tbtableSize) {
		if (_tbtable) {
			delete[] _tbtable;
		}
		_tbtableSize = bufferSize * 2; /*double the size*/
		_tbtable = new uint8_t[_tbtableSize];
		if (_tbtable == NULL) {
			Utils::exit("Memory allocation failed in function %s line %d\n",
					__FUNCTION__, __LINE__);
		}
	}

	/*build a query profile for s1*/
	bufferSize = 5 * s1Length;
	if (bufferSize > _profileSize) {
		if (_profile) {
			delete[] _profile;
		}
		_profileSize = bufferSize * 2;
		_profile = new int8_t[_profileSize];
		if (_profile == NULL) {
			Utils::exit("Memory allocation failed in function %s line %d\n",
					__FUNCTION__, __LINE__);
		}
	}
	profile = _profile;
	for (row = 0; row < 5; ++row) {
		pmatrix = _smat[row];
		for (col = 0; col < s1Length; ++col) {
			*profile++ = pmatrix[s1[col]];
		}
	}

	/*initialize the score buffer*/
	for (col = 0; col < s1Length; ++col) {
		_heBuffer[col] = make_int2(0, ALIGN_MIN_SCORE);
	}
	/*perform the semi-global alignment*/
	score = ALIGN_MIN_SCORE;
	cibest = cjbest = 0;
	for (row = 0; row < s2Length; ++row) {

		/*test if it is at the border*/
		border = row == s2Length - 1;

		/*initialilze the data*/
		h = p = 0;
		f = ALIGN_MIN_SCORE;

		/*get a row from the query profile*/
		pmatrix = _profile + s2[row] * s1Length;

		/*get the alignment backtracking matrix pointer*/
		pdata = _tbtable + row * _tbRowSize;
		for (col = 0; col < s1Length; ++col) {
			if ((h = h - _goe) > (f = f - _gext)) {
				f = h;
			}
			if ((h = _heBuffer[col].x - _goe)
					> (e = _heBuffer[col].y - _gext)) {
				e = h;
			}

			dir = ALIGN_DIR_DIAGONAL;
			h = p + pmatrix[col];
			//Utils::log("p: %d pmatrix %d\n", p, pmatrix[col]);
			if (h < f) {
				h = f;
				dir = ALIGN_DIR_LEFT;
			}
			if (h < e) {
				h = e;
				dir = ALIGN_DIR_UP;
			}

			p = _heBuffer[col].x;
			_heBuffer[col] = make_int2(h, e);

			/*check the last row and column*/
			if (border && col == s1Length - 1 && score < h) {
				score = h;
				cibest = row;
				cjbest = col;
			}

			/*save the alignment move*/
			setAlignMove(pdata, col, dir);
		}
	}
	//Utils::log("score %d cibest %d cjbest %d\n", score, cibest, cjbest);

	/*remove free end gaps*/
	int lastOp = ALIGN_DIR_STOP;
	do {
		/*check the coordinates*/
		dir = (cibest < 0 || cjbest < 0) ?
				ALIGN_DIR_STOP : getAlignMove(_tbtable + cibest * _tbRowSize, cjbest);

		if (dir == ALIGN_DIR_DIAGONAL || dir == ALIGN_DIR_STOP) {
			break;
		}
		/*check if the direction keeps the same*/
		if (lastOp == ALIGN_DIR_STOP) {
			lastOp = dir;
		}
		if (lastOp != dir) {
			break;
		}
		/*adjust the coordinates*/
		cibest -= dir == ALIGN_DIR_UP ? 1 : 0;
		cjbest -= dir == ALIGN_DIR_LEFT ? 1 : 0;
	} while (1);

	/*convert the alignment to CIGAR format*/
	CigarAlign* align = align2cigar(s1, s2, s1Length, s2Length, _tbtable, cibest, cjbest,
					score);
	if(align){
		aligns.push_back(align);
	}

#ifdef ALIGN_DEBUG
	traceback(s1, s2, s1Length, s2Length, _tbtable, cibest, cjbest, score);
#endif

	return aligns;
}
int32_t Aligner::traceback(uint8_t* s1, uint8_t* s2, int32_t s1Length,
		int32_t s2Length, uint8_t* tbtable, int32_t bandBufferSize,
		int32_t cibest, int32_t cjbest, int32_t alignScore, int32_t low) {
	bool done = false;
	int32_t dir;
	uint8_t c1, c2;
	int32_t ti, tj;
	int32_t row = cibest;
	;
	int32_t column = cjbest;
	int32_t alignLength = 0;
	int32_t numMismatches = 0;
	int32_t numIndels = 0;
	uint8_t *s1align, *s2align, *smalign;

	/*trace-back from the cell with the best score to get the banded optimal local CigarAlign*/
	s1align = new uint8_t[s1Length + bandBufferSize];
	if (s1align == NULL) {
		Utils::exit("Memory allocation failed in function %s line %d\n",
				__FUNCTION__, __LINE__);
	}
	s2align = new uint8_t[s2Length + bandBufferSize];
	if (s2align == NULL) {
		Utils::exit("Memory allocation failed in function %s line %d\n",
				__FUNCTION__, __LINE__);
	}

	smalign = new uint8_t[s2Length + bandBufferSize];
	if (smalign == NULL) {
		Utils::exit("Memory allocation failed in function %s line %d\n",
				__FUNCTION__, __LINE__);
	}
#ifdef ALIGN_DEBUG
	Utils::log("local alignment (s1 of length %d):\n", s1Length);
	for (int i = 0; i < s1Length; ++i) {
		Utils::log("%c", "ACGTN"[s1[i]]);
	}
	Utils::log("\n");

	Utils::log("local alignment (s2 of length %d):\n", s2Length);
	for (int i = 0; i < s2Length; ++i) {
		Utils::log("%c", "ACGTN"[s2[i]]);
	}
	Utils::log("\n");
	Utils::log("------------------------------------------\n");
#endif

	/*trace back to get the CigarAlign*/
	do {
		/*get the transformed coordinate in the CigarAlign matrix*/
		ti = row - 1;
		tj = (column + low + row - 1) - 1;

		/*check the direction of the alignment*/
		dir = tbtable[row * bandBufferSize + column];
		switch (dir) {
		case ALIGN_DIR_STOP:
			done = true;
			break;

		case ALIGN_DIR_DIAGONAL:
			/*save the aligned bases*/
			c1 = s1[ti];
			c2 = s2[tj];
			s1align[alignLength] = c1;
			s2align[alignLength] = c2;

			/*increase the number of mismatches*/
			if (c1 != c2) {
				++numMismatches;
				smalign[alignLength] = ' ';
			} else {
				smalign[alignLength] = '|';
			}
			++alignLength;

			/*adjusting the row and column*/
			--row;
			break;

		case ALIGN_DIR_UP: /*a deletion in s1*/
			c1 = s1[ti];
			s1align[alignLength] = c1;
			s2align[alignLength] = GAP_BASE;
			smalign[alignLength] = ' ';
			++alignLength;

			/*increase the number gaps*/
			++numIndels;

			/*adjust the row and column*/
			--row;
			++column;

			break;
		case ALIGN_DIR_LEFT: /*an insertion in s2*/
			c2 = s2[tj];
			s1align[alignLength] = GAP_BASE;
			s2align[alignLength] = c2;
			smalign[alignLength] = ' ';
			++alignLength;

			/*increase the number of gaps*/
			++numIndels;

			/*adjust the row and column*/
			--column;
			break;

		default:
			Utils::exit(
					"Unexpected value (%d) while tracing back the CigarAlign\n",
					dir);
			break;
		}
	} while (!done);

	/*NO local CigarAlign*/
	if (alignLength == 0) {
		delete[] smalign;
		delete[] s2align;
		delete[] s1align;
		return 0;
	}

	/*print out the CigarAlign*/
	int32_t numLines = (alignLength + 59) / 60;
	for (int32_t line = 0, index = alignLength - 1; line < numLines;
			++line, index -= 60) {
		for (int32_t i = 0; i < 60 && index - i >= 0; ++i) {
			fputc(decode(s1align[index - i]), stderr);
		}
		fputc('\n', stderr);
		for (int32_t i = 0; i < 60 && index - i >= 0; ++i)
		{
			fputc(smalign[index - i], stderr);
		}
		fputc('\n', stderr);
		for (int32_t i = 0; i < 60 && index - i >= 0; ++i)
		{
			fputc(decode(s2align[index - i]), stderr);
		}
		fputc('\n', stderr);
	}

	delete[] smalign;
	delete[] s2align;
	delete[] s1align;

	Utils::log("#mismatches: %d #indels: %d #alignScore: %d #alignLength: %d\n",
			numMismatches, numIndels, alignScore, alignLength);
	return numMismatches + numIndels;
}

int32_t Aligner::traceback(uint8_t* s1, uint8_t* s2, int32_t s1Length,
		int32_t s2Length, uint8_t* tbtable, int32_t cibest, int32_t cjbest,
		int32_t alignScore) {
	int32_t dir;
	uint8_t c1, c2;
	int32_t row = cibest;
	int32_t column = cjbest;
	int32_t alignLength = 0;
	int32_t numMismatches = 0;
	int32_t numIndels = 0;
	uint8_t *s1align, *s2align, *smalign;

	/*trace-back from the cell with the best score to get the banded optimal local CigarAlign*/
	s1align = new uint8_t[s1Length + s2Length];
	if (s1align == NULL) {
		Utils::exit("Memory allocation failed in function %s line %d\n",
				__FUNCTION__, __LINE__);
	}
	s2align = new uint8_t[s1Length + s2Length];
	if (s2align == NULL) {
		Utils::exit("Memory allocation failed in function %s line %d\n",
				__FUNCTION__, __LINE__);
	}

	smalign = new uint8_t[s1Length + s2Length];
	if (smalign == NULL) {
		Utils::exit("Memory allocation failed in function %s line %d\n",
				__FUNCTION__, __LINE__);
	}

#ifdef ALIGN_DEBUG
	Utils::log("local alignment (s1 of length %d):\n", s1Length);
	for (int i = 0; i < s1Length; ++i) {
		Utils::log("%c", "ACGTN"[s1[i]]);
	}
	Utils::log("\n");

	Utils::log("local alignment (s2 of length %d):\n", s2Length);
	for (int i = 0; i < s2Length; ++i) {
		Utils::log("%c", "ACGTN"[s2[i]]);
	}
	Utils::log("\n");
	Utils::log("------------------------------------------\n");
#endif

	/*trace back to get the CigarAlign*/
	do {
		if (row == -1 || column == -1) {
			break;
		}
		dir = getAlignMove(tbtable + row * _tbRowSize, column);
		switch (dir) {
		case ALIGN_DIR_DIAGONAL:
			/*save the aligned bases*/
			c1 = s1[column];
			c2 = s2[row];
			s1align[alignLength] = c1;
			s2align[alignLength] = c2;

			/*increase the number of mismatches*/
			if (c1 != c2) {
				++numMismatches;
				smalign[alignLength] = ' ';
			} else {
				smalign[alignLength] = '|';
			}
			++alignLength;
			--row;
			--column;
			break;

		case ALIGN_DIR_UP: /*a deletion in s2*/
			c2 = s2[row];
			s1align[alignLength] = GAP_BASE;
			s2align[alignLength] = c2;
			smalign[alignLength] = ' ';
			++alignLength;

			/*increase the number gaps*/
			++numIndels;
			--row;
			break;
		case ALIGN_DIR_LEFT: /*an insertion in s1*/
			c1 = s1[column];
			s1align[alignLength] = c1;
			s2align[alignLength] = GAP_BASE;
			smalign[alignLength] = ' ';
			++alignLength;

			/*increase the number of gaps*/
			++numIndels;
			--column;
			break;

		default:
			Utils::exit(
					"Unexpected value (%d) while tracing back the CigarAlign\n",
					dir);
			break;
		}
	} while (1);

	/*NO local CigarAlign*/
	if (alignLength == 0) {
		delete[] smalign;
		delete[] s2align;
		delete[] s1align;
		return 0;
	}

	/*print out the CigarAlign*/
	int32_t numLines = (alignLength + 59) / 60;
	for (int32_t line = 0, index = alignLength - 1; line < numLines;
			++line, index -= 60) {
		for (int32_t i = 0; i < 60 && index - i >= 0; ++i) {
			fputc(decode(s1align[index - i]), stderr);
		}
		fputc('\n', stderr);
		for (int32_t i = 0; i < 60 && index - i >= 0; ++i)
		{
			fputc(smalign[index - i], stderr);
		}
		fputc('\n', stderr);
		for (int32_t i = 0; i < 60 && index - i >= 0; ++i)
		{
			fputc(decode(s2align[index - i]), stderr);
		}
		fputc('\n', stderr);
	}

	delete[] smalign;
	delete[] s2align;
	delete[] s1align;

	Utils::log("#mismatches: %d #indels: %d #alignScore: %d #alignLength: %d\n",
			numMismatches, numIndels, alignScore, alignLength);
	return numMismatches + numIndels;
}
CigarAlign* Aligner::align2cigar(uint8_t* s1, uint8_t* s2, int32_t s1Length,
		int32_t s2Length, uint8_t* tbtable, int32_t bandBufferSize,
		int32_t cibest, int32_t cjbest, int32_t alignScore, int32_t low,
		bool rowwiseQuery) {
	bool done = false;
	int32_t dir;
	uint8_t c1, c2;
	int32_t ti, tj;
	int32_t row = cibest, prow = cibest;
	int32_t column = cjbest, pcolumn = cjbest;
	int32_t alignLength = 0;
	int32_t numMismatches = 0, numIndels = 0, numGaps = 0;
	uint8_t op = ' ';
	uint8_t lastOp = ALIGN_DIR_STOP;
	uint32_t numOps = 0;
	CigarAlign* align;
	uint32_t* cigar;
	int32_t cigarNum, cigarSize;
	int32_t start1, end1, start2, end2;

	/*allocate space for cigar for s1, assuming that each entry takes 8 bytes*/
	cigarNum = 0;
	cigarSize = 128;
	cigar = new uint32_t[cigarSize];
	if (cigar == NULL) {
		Utils::exit("Memory allocation failed in function %s line %d\n",
				__FUNCTION__, __LINE__);
	}

	/*trace back to get the CigarAlign*/
	while (!done) {
		/*get the transformed coordinate in the CigarAlign matrix*/
		ti = row - 1;
		tj = (column + low + row - 1) - 1;

		/*check the direction of the alignment*/
		dir = tbtable[row * bandBufferSize + column];
		switch (dir) {
		case ALIGN_DIR_STOP:
			done = true;
			/*record the last operations*/
			if (numOps > 0) {
				/*resize the memory*/
				if (cigarNum >= cigarSize) {
					cigarSize += 128;
					uint32_t* buffer = new uint32_t[cigarSize];
					if (buffer == NULL) {
						Utils::exit(
								"Memory allocation failed in function %s line %d\n",
								__FUNCTION__, __LINE__);
					}
					memcpy(buffer, cigar, cigarNum * sizeof(cigar[0]));
					delete[] cigar;
					cigar = buffer;
				}
				cigar[cigarNum++] = (numOps << 2) | op;
				numGaps += op != ALIGN_OP_M ? 1 : 0;
			}
			break;

		case ALIGN_DIR_DIAGONAL:
			/*save the aligned bases*/
			c1 = s1[ti];
			c2 = s2[tj];
			++alignLength;

			/*increase the number of mismatches*/
			if (c1 != c2) {
				++numMismatches;
			}
			/*check the operation type*/
			if (lastOp == ALIGN_DIR_DIAGONAL) {
				++numOps;
			} else {
				/*save the previous operation*/
				if (numOps > 0) {
					/*resize the memory*/
					if (cigarNum >= cigarSize) {
						cigarSize += 1024;
						uint32_t* buffer = new uint32_t[cigarSize];
						if (buffer == NULL) {
							Utils::exit(
									"Memory allocation failed in function %s line %d\n",
									__FUNCTION__, __LINE__);
						}
						memcpy(buffer, cigar, cigarNum * sizeof(cigar[0]));
						delete[] cigar;
						cigar = buffer;
					}
					cigar[cigarNum++] = (numOps << 2) | op;
					numGaps += op != ALIGN_OP_M ? 1 : 0;
				}

				/*change to a new operation*/
				numOps = 1;
				op = ALIGN_OP_M;
				lastOp = ALIGN_DIR_DIAGONAL;
			}

			/*adjusting the row and column*/
			prow = row;
			pcolumn = column;
			--row;
			break;

		case ALIGN_DIR_UP: /*an insertion in s1*/
			++alignLength;

			/*increase the number gaps*/
			++numIndels;

			/*check the operation type*/
			if (lastOp == ALIGN_DIR_UP) {
				++numOps;
			} else {
				/*save the previous operation*/
				if (numOps > 0) {
					/*resize the memory*/
					if (cigarNum >= cigarSize) {
						cigarSize += 1024;
						uint32_t* buffer = new uint32_t[cigarSize];
						if (buffer == NULL) {
							Utils::exit(
									"Memory allocation failed in function %s line %d\n",
									__FUNCTION__, __LINE__);
						}
						memcpy(buffer, cigar, cigarNum * sizeof(cigar[0]));
						delete[] cigar;
						cigar = buffer;
					}
					cigar[cigarNum++] = (numOps << 2) | op;
					numGaps += op != ALIGN_OP_M ? 1 : 0;
				}

				/*change to a new operation*/
				numOps = 1;
				op = ALIGN_OP_I;
				lastOp = ALIGN_DIR_UP;
			}

			/*adjust the row and column*/
			prow = row;
			pcolumn = column;
			--row;
			++column;

			break;
		case ALIGN_DIR_LEFT: /*a deletion in s1*/
			++alignLength;

			/*increase the number of gaps*/
			++numIndels;

			/*check the operation type*/
			if (lastOp == ALIGN_DIR_LEFT) {
				++numOps;
			} else {
				/*save the previous operation*/
				if (numOps > 0) {
					/*resize the memory*/
					if (cigarNum >= cigarSize) {
						cigarSize += 1024;
						uint32_t* buffer = new uint32_t[cigarSize];
						if (buffer == NULL) {
							Utils::exit(
									"Memory allocation failed in function %s line %d\n",
									__FUNCTION__, __LINE__);
						}
						memcpy(buffer, cigar, cigarNum * sizeof(cigar[0]));
						delete[] cigar;
						cigar = buffer;
					}
					cigar[cigarNum++] = (numOps << 2) | op;
					numGaps += op != ALIGN_OP_M ? 1 : 0;
				}

				/*change to a new operation*/
				numOps = 1;
				op = ALIGN_OP_D;
				lastOp = ALIGN_DIR_LEFT;
			}

			/*adjust the row and column*/
			prow = row;
			pcolumn = column;
			--column;
			break;

		default:
			Utils::exit(
					"Unexpected value (%d) while tracing back the CigarAlign\n",
					dir);
			break;
		}
	}

	/*NO local CigarAlign*/
	if (alignLength == 0) {
		delete[] cigar;
		return NULL;
	}

	/*create an cigar-based CigarAlign object for s1*/
	start1 = prow - 1;
	start2 = (pcolumn + low + prow - 1) - 1;
	end1 = cibest - 1;
	end2 = (cjbest + low + cibest - 1) - 1;

	if (rowwiseQuery == true) {
		/*create alignment for the first sequence*/
		align = new CigarAlign(REALIGN_TYPE_LOCAL, alignLength, start1, end1,
				start2, end2, numMismatches, numIndels, numGaps, cigar,
				cigarNum, alignScore, rowwiseQuery);
		//Utils::log("rowwiseQuery %d start1 %d end1 %d start2 %d end2 %d\n", rowwiseQuery, start1, end1, start2, end2);
	} else {
		/*create alignment for the second sequence*/
		align = new CigarAlign(REALIGN_TYPE_LOCAL, alignLength, start2, end2,
				start1, end1, numMismatches, numIndels, numGaps, cigar,
				cigarNum, alignScore, rowwiseQuery);
		//Utils::log("rowwiseQuery %d start1 %d end1 %d start2 %d end2 %d\n", rowwiseQuery, start2, end2, start1, end1);
	}

	/*release buffer*/
	delete[] cigar;

	return align;
}

/*trace-back of semi-global alignment*/
CigarAlign* Aligner::align2cigar(uint8_t* s1, uint8_t* s2, int32_t s1Length,
		int32_t s2Length, uint8_t* tbtable, int32_t cibest, int32_t cjbest,
		int32_t alignScore) {
	int32_t dir;
	uint8_t c1, c2;
	int32_t row = cibest, prow = cibest;
	int32_t column = cjbest, pcolumn = cjbest;
	int32_t alignLength = 0;
	int32_t numMismatches = 0, numIndels = 0, numGaps = 0;
	uint8_t op = ALIGN_OP_S;
	uint8_t lastOp = ALIGN_DIR_STOP;
	uint32_t numOps = 0;
	CigarAlign* align;
	uint32_t* cigar;
	int32_t cigarNum, cigarSize;

	/*recaluclate alignment score*/
	alignScore = 0;

	/*allocate space for cigar for s1, assuming that each entry takes 8 bytes*/
	cigarNum = 0;
	cigarSize = 128;
	cigar = new uint32_t[cigarSize];
	if (cigar == NULL) {
		Utils::exit("Memory allocation failed in function %s line %d\n",
				__FUNCTION__, __LINE__);
	}

	/*trace back to get the CigarAlign*/
	bool isDiagnal;
	do {
		/*check the coordinates*/
		dir = (row < 0 || column < 0) ?
				ALIGN_DIR_STOP : getAlignMove(tbtable + row * _tbRowSize, column);
		isDiagnal = dir == ALIGN_DIR_DIAGONAL;

		/*count and save the number of operations*/
		if (lastOp == dir) {
			++numOps;
		} else {
			if (numOps) {
				/*resize the memory*/
				if (cigarNum >= cigarSize) {
					cigarSize += 128;
					uint32_t* buffer = new uint32_t[cigarSize];
					if (buffer == NULL) {
						Utils::exit(
								"Memory allocation failed in function %s line %d\n",
								__FUNCTION__, __LINE__);
					}
					memcpy(buffer, cigar, cigarNum * sizeof(cigar[0]));
					delete[] cigar;
					cigar = buffer;
				}
				cigar[cigarNum++] = (numOps << 2) | op;
				numGaps += op != ALIGN_OP_M ? 1 : 0;
			}
			/*reset the values*/
			numOps = 1;
			lastOp = dir;
			op = isDiagnal ?
					ALIGN_OP_M :
					(dir == ALIGN_DIR_LEFT ? ALIGN_OP_I : ALIGN_OP_D);
		}

		/*exit the backtracking*/
		if (dir == ALIGN_DIR_STOP) {
			break;
		}

		/*calculate the edit distances*/
		if (isDiagnal) {
			c1 = s1[column];
			c2 = s2[row];
		}
		numMismatches += isDiagnal ? c1 != c2 : 0;
		numIndels += isDiagnal ? 0 : 1;
		if(isDiagnal){
			alignScore += c1 == c2 ? _match : _mismatch;
		}else{
			alignScore -= numOps == 1 ? _goe : _gext;
		}
		++alignLength;

		/*adjust the coordinates*/
		prow = row;
		pcolumn = column;
		column -= (isDiagnal || dir == ALIGN_DIR_LEFT);
		row -= (isDiagnal || dir == ALIGN_DIR_UP);
	} while (1);

	/*NO local CigarAlign*/
	if (alignLength == 0) {
		delete[] cigar;
		return NULL;
	}
	/*create alignment for the second sequence*/
	align = new CigarAlign(REALIGN_TYPE_SEMIGLOBAL, alignLength, pcolumn, cjbest,
			prow, cibest, numMismatches, numIndels, numGaps, cigar, cigarNum,
			alignScore, true);

	/*release buffer*/
	delete[] cigar;

	return align;
}

/***************************************************************************
 *SSE-accelerated Smith-Waterman optimal local alignment score computation
 **************************************************************************/

int32_t Aligner::search7(BYTE** qtable, BYTE gap_open_penalty,
		BYTE gap_extend_penalty, BYTE * score_matrix, BYTE * dprofile,
		BYTE * hearray, int32_t qlen, int32_t numSeqs, int32_t *seqOffsets,
		uint8_t *sequences, AlignScore* scores) {

	int32_t maxScore = 0;
	__m128i S, Q, R, T, M, Z, T0;
	__m128i *hep, **qp;
	BYTE * d_begin[N16_CHANNELS];

	__m128i dseqalloc[CDEPTH];

	BYTE * dseq = (BYTE*) &dseqalloc;
	BYTE zero;

	long seq_id[N16_CHANNELS];
	long next_id = 0;
	int32_t done;

	memset(hearray, 0x80, qlen * 32);

	Z = _mm_set_epi8(0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80,
			0x80, 0x80, 0x80, 0x80, 0x80, 0x80);
	T0 = _mm_set_epi8(0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
			0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x80);
	Q = _mm_set_epi8(gap_open_penalty, gap_open_penalty, gap_open_penalty,
			gap_open_penalty, gap_open_penalty, gap_open_penalty,
			gap_open_penalty, gap_open_penalty, gap_open_penalty,
			gap_open_penalty, gap_open_penalty, gap_open_penalty,
			gap_open_penalty, gap_open_penalty, gap_open_penalty,
			gap_open_penalty);
	R = _mm_set_epi8(gap_extend_penalty, gap_extend_penalty, gap_extend_penalty,
			gap_extend_penalty, gap_extend_penalty, gap_extend_penalty,
			gap_extend_penalty, gap_extend_penalty, gap_extend_penalty,
			gap_extend_penalty, gap_extend_penalty, gap_extend_penalty,
			gap_extend_penalty, gap_extend_penalty, gap_extend_penalty,
			gap_extend_penalty);
	zero = 0;
	done = 0;

	S = Z;

	hep = (__m128i *) hearray;
	qp = (__m128i **) qtable;

	for (int c = 0; c < N16_CHANNELS; c++) {
		d_begin[c] = &zero;
		seq_id[c] = -1;
	}

	int easy = 0;

	while (1) {
		if (easy) {
			// fill all channels

			for (int c = 0; c < N16_CHANNELS; c++) {
				for (int j = 0; j < CDEPTH; j++) {
					BYTE v = *(d_begin[c]);
					dseq[N16_CHANNELS * j + c] = v;
					if (v)
						d_begin[c]++;
				}
				if (!*(d_begin[c]))
					easy = 0;
			}

			if (_haveSSSE3) {
				dprofile_shuffle7(dprofile, score_matrix, dseq);
			} else {
				dprofile_fill7(dprofile, score_matrix, dseq);
			}

			donormal7(&S, hep, qp, &Q, &R, qlen, &Z);
		} else {
			// One or more sequences ended in the previous block
			// We have to switch over to a new sequence

			easy = 1;

			M = _mm_setzero_si128();
			T = T0;
			for (int c = 0; c < N16_CHANNELS; c++) {
				if (*(d_begin[c])) {
					// this channel has more sequence

					for (int j = 0; j < CDEPTH; j++) {
						BYTE v = *(d_begin[c]);
						dseq[N16_CHANNELS * j + c] = v;
						if (v)
							d_begin[c]++;
					}
					if (!*(d_begin[c]))
						easy = 0;
				} else {
					// sequence in channel c ended
					// change of sequence

					M = _mm_xor_si128(M, T);

					long cand_id = seq_id[c];

					if (cand_id >= 0) {
						// save score
						long score = ((BYTE*) &S)[c] - 0x80;
						scores[cand_id]._score = score;
						if (score > maxScore) {
							maxScore = score;
						}
						done++;
					}
					if (next_id < numSeqs) {
						// get next sequence
						seq_id[c] = next_id;
						d_begin[c] = sequences + seqOffsets[next_id];
						next_id++;

						// fill channel
						for (int j = 0; j < CDEPTH; j++) {
							BYTE v = *(d_begin[c]);
							dseq[N16_CHANNELS * j + c] = v;
							if (v)
								d_begin[c]++;
						}
						if (!*(d_begin[c]))
							easy = 0;
					} else {
						// no more sequences, empty channel
						seq_id[c] = -1;
						d_begin[c] = &zero;
						for (int j = 0; j < CDEPTH; j++)
							dseq[N16_CHANNELS * j + c] = 0;
					}

				}

				T = _mm_slli_si128(T, 1);
			}

			if (done == numSeqs)
				break;

			if (_haveSSSE3) {
				dprofile_shuffle7(dprofile, score_matrix, dseq);
			} else {
				dprofile_fill7(dprofile, score_matrix, dseq);
			}

			domasked7(&S, hep, qp, &Q, &R, qlen, &Z, &M);
		}
	}
	return maxScore;
}
void Aligner::search16(WORD** qtable, WORD gap_open_penalty,
		WORD gap_extend_penalty, WORD * score_matrix, WORD * dprofile,
		WORD * hearray, int32_t qlen, int32_t numSeqs, int32_t *seqOffsets,
		uint8_t *sequences, AlignScore* scores) {

	__m128i S, Q, R, T, M, Z, T0;
	__m128i *hep, **qp;
	BYTE * d_begin[N8_CHANNELS];

	__m128i dseqalloc[CDEPTH];

	BYTE * dseq = (BYTE *) &dseqalloc;
	BYTE zero;

	int32_t seq_id[N8_CHANNELS];
	int32_t next_id = 0;
	int32_t done;

	Z = _mm_set_epi16(0x8000, 0x8000, 0x8000, 0x8000, 0x8000, 0x8000, 0x8000,
			0x8000);
	T0 = _mm_set_epi16(0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000,
			0x8000);
	Q = _mm_set_epi16(gap_open_penalty, gap_open_penalty, gap_open_penalty,
			gap_open_penalty, gap_open_penalty, gap_open_penalty,
			gap_open_penalty, gap_open_penalty);
	R = _mm_set_epi16(gap_extend_penalty, gap_extend_penalty,
			gap_extend_penalty, gap_extend_penalty, gap_extend_penalty,
			gap_extend_penalty, gap_extend_penalty, gap_extend_penalty);

	zero = 0;
	done = 0;

	S = Z;

	hep = (__m128i *) hearray;
	qp = (__m128i **) qtable;

	for (int a = 0; a < qlen; a++) {
		hep[2 * a] = Z;
		hep[2 * a + 1] = Z;
	}

	for (int c = 0; c < N8_CHANNELS; c++) {
		d_begin[c] = &zero;
		seq_id[c] = -1;
	}

	int easy = 0;

	while (1) {
		if (easy) {
			for (int c = 0; c < N8_CHANNELS; c++) {
				for (int j = 0; j < CDEPTH; j++) {
					BYTE v = *(d_begin[c]);
					dseq[N8_CHANNELS * j + c] = v;
					if (v)
						d_begin[c]++;
				}
				if (!*(d_begin[c]))
					easy = 0;
			}

			dprofile_fill16(dprofile, score_matrix, dseq);

			donormal16(&S, hep, qp, &Q, &R, qlen, &Z);

		} else {

			easy = 1;

			M = _mm_setzero_si128();
			T = T0;

			for (int c = 0; c < N8_CHANNELS; c++) {
				if (*(d_begin[c])) {
					for (int j = 0; j < CDEPTH; j++) {
						BYTE v = *(d_begin[c]);
						dseq[N8_CHANNELS * j + c] = v;
						if (v)
							d_begin[c]++;
					}

					if (!*(d_begin[c]))
						easy = 0;

				} else {
					M = _mm_xor_si128(M, T);

					long cand_id = seq_id[c];

					if (cand_id >= 0) {
						int score = ((WORD*) &S)[c] - 0x8000;
						/*save the alignment score*/
						scores[cand_id]._score = score;
						done++;
					}
#ifndef ONLY_8CHANNEL
					/*find the next non-processed sequence*/
					for (;
							next_id < numSeqs
									&& scores[next_id]._score < _scorelimit7;
							++next_id) {
						done++;
					}
#endif
					if (next_id < numSeqs) {
						seq_id[c] = next_id;
						d_begin[c] = sequences + seqOffsets[next_id];
						next_id++;

						for (int j = 0; j < CDEPTH; j++) {
							BYTE v = *(d_begin[c]);
							dseq[N8_CHANNELS * j + c] = v;
							if (v)
								d_begin[c]++;
						}
						if (!*(d_begin[c]))
							easy = 0;
					} else {
						seq_id[c] = -1;
						d_begin[c] = &zero;
						for (int j = 0; j < CDEPTH; j++)
							dseq[N8_CHANNELS * j + c] = 0;
					}
				}
				T = _mm_slli_si128(T, 2);
			}

			if (done == numSeqs)
				break;

			dprofile_fill16(dprofile, score_matrix, dseq);

			domasked16(&S, hep, qp, &Q, &R, qlen, &Z, &M);
		}
	}
}
