/*
 * Mapping.h
 *
 *  Created on: Jan 18, 2012
 *      Author: yongchao
 */

#ifndef MAPPING_H_
#define MAPPING_H_
#include "Macros.h"

struct Mapping
{
	Mapping() {
		_align = NULL;
		_data = NULL;
		_seqLength = 0;
	}
	Mapping(CigarAlign* align, int64_t position, int64_t gposition, int strand,
			int genomeIndex, int mapQual, int mapQualBase) {
		_align = NULL;
		_data = NULL;
		_seqLength = 0;

		init(align, position, gposition, strand, genomeIndex, mapQual,
				mapQualBase);
	}
	~Mapping() {
		if (_align) {
			delete _align;
		}
		if (_data) {
			delete[] _data;
		}
	}
	inline void setData(uint8_t* data, uint32_t seqLength) {
		_data = data;
		_seqLength = seqLength;
	}
	inline void copyData(Mapping* mapping) {
		if (_data) {
			delete[] _data;
		}
		_seqLength = mapping->_seqLength;
		_data = new uint8_t[_seqLength];
		if (!_data) {
			Utils::exit("Memory allocation (%u) failed at line %d in file %s\n",
					_seqLength, __LINE__, __FILE__);
		}
		if (mapping->_data) {
			memcpy(_data, mapping->_data, _seqLength * 3);
		} else {
			Utils::exit("There MUST be something wrong with the program\n");
		}
	}
	inline void init(CigarAlign* align, int64_t position, int64_t gposition,
			int strand, int genomeIndex, int mapQual, int mapQualBase) {
		if (_align) {
			delete _align;
		}
		_align = align;
		_position = position;
		_gposition = gposition;
		_genomeIndex = genomeIndex;
		_strand = strand;
		_mapQual = mapQual;
		_mapQualBase = mapQualBase;
	}

	/*member variables*/
	CigarAlign* _align;
	int64_t _position; /*mapping position on the chromosome*/
	int64_t _gposition; /*mapping positon on the whole genome (concatenated chromosomes*/
	int _genomeIndex;
	uint32_t _strand :1;
	uint32_t _mapQualBase :16;
	uint32_t _mapQual :15;

	/*for converted color-space read*/
	uint8_t* _data;
	uint32_t _seqLength; /*sequence length*/
};

struct MappingPair
{
	MappingPair() {
		_score = -1;
	}
	MappingPair(Mapping* left, Mapping* right, float score) {
		_left = left;
		_right = right;
		_score = score;
		_editDist = _left->_align->getEditDistance()
				+ _right->_align->getEditDistance();
	}

	Mapping *_left;
	Mapping *_right;
	float _score;
	int32_t _editDist;
};

class MappingPairComp {
public:
	inline bool operator()(MappingPair& a, MappingPair &b) const {
		if (a._score < b._score
				|| (a._score == b._score && a._editDist > b._editDist)) {
			return true;
		}
		return false;
	}
};

#endif /* MAPPING_H_ */
