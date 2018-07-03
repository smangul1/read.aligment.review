/*
 *  macros.h
 *
 *  Created on: Dec 23, 2011
 *      Author: yongchao
 */

#ifndef MACROS_H_
#define MACROS_H_
#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <stdint.h>
#include <assert.h>
#include <limits.h>
#include <zlib.h>
#include <math.h>
#include <pthread.h>
#include <iterator>
#include <algorithm>
using namespace std;

/*the version and name*/
#define PROGRAM_NAME	"cushaw3"
#define PROGRAM_VERSION "3.0.3"

/*for SSE2 implementations*/
#ifdef HAVE_SSSE3
#include <tmmintrin.h>
#else
#include <emmintrin.h>
#endif

typedef unsigned int UINT32;
typedef unsigned short WORD;
typedef unsigned char BYTE;
typedef BYTE VECTOR[16];

//#define max(a, b)	(((a) < (b)) ? (b) : (a))
//#define min(a, b)	(((a) < (b)) ? (a) : (b))

/*maximal number of top alignments*/
#define MAX_MULTI_ALIGN			1000

/*macros for the alignment path trace-back*/
#define ALIGN_DIR_STOP				0
#define	ALIGN_DIR_DIAGONAL			1
#define ALIGN_DIR_UP				2
#define	ALIGN_DIR_LEFT				3
#define ALIGN_MIN_SCORE				-999999

//SAM format
#define SAM_FPD   1 // paired
#define SAM_FPP   2 // properly paired
#define SAM_FSU   4 // self-unmapped
#define SAM_FMU   8 // mate-unmapped
#define SAM_FSR  16 // self on the reverse strand
#define SAM_FMR  32 // mate on the reverse strand
#define SAM_FR1  64 // this is read one
#define SAM_FR2 128 // this is read two
#define SAM_FSC 256 // secondary alignment
#define	DEFAULT_MAX_MAP_QUAL	250
#define SW_MAP_QUALITY_SCORE	0.5

//single and paired end alignment
#define MAX_SEQ_LENGTH					4096		/*for higher values, pre-computed tables are not used*/
#define GLOBAL_MIN_SEED_SIZE			8
#define GLOBAL_MAX_SEED_SIZE			49
#define GLOBAL_MAX_NUM_SEED_REPEATS		2000
#define GLOBAL_MIN_NUM_SEED_PAIRS		1000000		/*1 million*/
#define GLOBAL_MAX_NUM_SEED_PAIRS		10000000	/*10 million*/

#define DEFAULT_SEQ_LENGTH_RATIO		0.80
#define DEFAULT_MIN_IDENTITY			0.90
#define DEFAULT_GLOBAL_MIN_IDENTITY		0.65

#define FILE_TYPE_FASTX   0
#define FILE_TYPE_BAM     1
#define FILE_TYPE_SAM     2
#define FILE_FORMAT_FASTQ   3
#define FILE_FORMAT_FASTA   4
#define FILE_FORMAT_BSAM    5

#define INS_SIZE_EST_MULTIPLE		0x10000

/*alignment types*/
#define HAVE_LOCAL_ALIGN				0
#define HAVE_LOCAL_SEMI_GLOBAL_ALIGN	1


/*pairing mode*/
#define PAIRED_END_READS		0
#define MATE_PAIRED_READS		1

/*alignment algorithm*/
#define REALIGN_TYPE_LOCAL			0
#define REALIGN_TYPE_SEMIGLOBAL		1

/*seed type*/
#define MAXIMAL_EXACT_MATCH_SEED	0
#define FIXED_LENGTH_STRIPED_SEED	1

#include "Structs.h"
#include "BwtMacros.h"
#endif /* MACROS_H_ */
