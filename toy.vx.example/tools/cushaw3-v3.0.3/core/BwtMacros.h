/*
 * BwtMacros.h
 *
 *  Created on: Dec 23, 2011
 *      Author: yongchao
 */

#ifndef BWTMACROS_H_
#define BWTMACROS_H_

#ifndef _BWT_DEFS_H
#define _BWT_DEFS_H

//Macros for burrows-wheeler transform
#define BWT_NUM_NUCLEOTIDE			4
#define BWT_NUM_OCC					(BWT_NUM_NUCLEOTIDE + 1)
#define BWT_OCC_INTERVAL			128		//must be multiples of 16
#define BWT_OCC_INTERVAL_SHIFT		7
#define BWT_OCC_INTERVAL_MASK		127
#define BWT_OCC_PTR_OFFSET			(BWT_NUM_NUCLEOTIDE + (BWT_OCC_INTERVAL >> 4))	//a 32-bit integer stores 16 bases
//for alignment information
#define BWT_MAX_SCORE				0x0FFFFE
#define BWT_SCORE_MASK				0x0FFFFF		//20 bits
#define BWT_STRAND_BIT				0x80000000
#define BWT_STRAND_BIT_MASK			0x80000000
#define BWT_STRAND_BIT_SHIFT		31			//one bit
#define BWT_SEQ_STRAND_BIT			0x40000000
#define BWT_SEQ_STRAND_BIT_SHIFT	30			//one bit
#define BWT_NUM_MISMATCH_SHIFT		20
#define BWT_NUM_MISMATCH_MASK		0x3ff		//10 bits
//sequence base encoding
#define ADENINE			0
#define CYTOSINE		1
#define GUANINE			2
#define THYMINE			3
#define	UNKNOWN_BASE	BWT_NUM_NUCLEOTIDE
#define GAP_BASE		(UNKNOWN_BASE + 1)

/*decode a base*/
#define decode(base)		("ACGTN-"[base])

#endif

#endif /* BWTMACROS_H_ */
