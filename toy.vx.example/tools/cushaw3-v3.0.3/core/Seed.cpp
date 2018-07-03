/*
 * Seed.cpp
 *
 *  Created on: Feb 3, 2012
 *      Author: yongchao
 */

#include "Seed.h"
#include "Utils.h"
bool operator>(const AlignScore& p, const AlignScore& q) {
	return (p._score < q._score
			|| (p._score == q._score && p._seedLength < q._seedLength));
}
bool operator<(const AlignScore& p, const AlignScore& q) {
	return (p._score > q._score
			|| (p._score == q._score && p._seedLength > q._seedLength));
}
bool operator==(const AlignScore& p, const AlignScore& q) {
	return (p._score == q._score && p._seedLength == q._seedLength);
}

bool operator>(const Seed& p, const Seed& q) {
	return (p._alignScore < q._alignScore
			|| (p._alignScore == q._alignScore && p._seedLength < q._seedLength));
}
bool operator<(const Seed& p, const Seed& q) {
	return (p._alignScore > q._alignScore
			|| (p._alignScore == q._alignScore && p._seedLength > q._seedLength));
}
bool operator==(const Seed& p, const Seed& q) {
	return (p._alignScore == q._alignScore && p._seedLength == q._seedLength);
}
