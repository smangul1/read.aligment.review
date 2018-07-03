/*
 * Seed.h
 *
 *  Created on: Jan 26, 2012
 *      Author: yongchao
 */

#ifndef SEED_H_
#define SEED_H_
#include "Macros.h"

struct Seed
{
	uint32_t _targetPosition;
	uint32_t _strand :1;
	uint32_t _alignScore : 31;
	uint32_t _queryPosition :16;
	uint32_t _seedLength :16;
	int32_t _genomeIndex;
};

struct SeedPair
{
	SeedPair()
	{
		_score = -1;
	}
	SeedPair(Seed& left, Seed& right, float score){
		_left = left;
		_right = right;
		_score = score;
	}
	Seed _left;
	Seed _right;
	float _score;
};

class SeedPairComp {
public:
  inline bool operator()(SeedPair& a, SeedPair &b) const {
    return (a._score < b._score) ? true : false;
  }
};

struct AlignScore{
	int32_t _score;
	uint32_t _seedLength :16;
	uint32_t _seedPosition :16;
	int32_t _seedIndex;
	uint32_t _targetPosition;
};

bool operator>(const AlignScore& p, const AlignScore& q);
bool operator<(const AlignScore& p, const AlignScore& q);
bool operator==(const AlignScore& p, const AlignScore& q);

bool operator>(const Seed& p, const Seed& q);
bool operator<(const Seed& p, const Seed& q);
bool operator==(const Seed& p, const Seed& q);
#endif /* SEED_H_ */
