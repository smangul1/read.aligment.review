/*
 * Structs.h
 *
 *  Created on: Dec 27, 2011
 *      Author: yongchao
 */

#ifndef STRUCTS_H_
#define STRUCTS_H_

#ifndef __CUDACC__

/*structure for int2*/
typedef struct {
	int32_t x;
	int32_t y;
} int2;
#define make_int2(x, y) (int2){x, y}

/*structure for uint2*/
typedef struct {
	uint32_t x;
	uint32_t y;
} uint2;
#define make_uint2(x, y) (uint2){x, y}

/*structure for uint3*/
typedef struct {
	uint32_t x;
	uint32_t y;
	uint32_t z;
} uint3;
#define make_uint3(x, y, z) (uint3){x, y, z}

/*structure for uint4*/
typedef struct {
	uint32_t x;
	uint32_t y;
	uint32_t z;
	uint32_t w;
} uint4;
#define make_uint4(x, y, z, w) (uint4){x, y, z, w}

#include "Seed.h"
#endif

#endif /* STRUCTS_H_ */
