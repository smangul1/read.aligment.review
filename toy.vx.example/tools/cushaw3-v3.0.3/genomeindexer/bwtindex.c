/* The MIT License

 Copyright (c) 2008 Genome Research Ltd (GRL).

 Permission is hereby granted, free of charge, to any person obtaining
 a copy of this software and associated documentation files (the
 "Software"), to deal in the Software without restriction, including
 without limitation the rights to use, copy, modify, merge, publish,
 distribute, sublicense, and/or sell copies of the Software, and to
 permit persons to whom the Software is furnished to do so, subject to
 the following conditions:

 The above copyright notice and this permission notice shall be
 included in all copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
 BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
 ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 SOFTWARE.
 */

/* Contact: Heng Li <lh3@sanger.ac.uk> */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <zlib.h>
#include <math.h>
#include "bntseq.h"
#include "bwt.h"
#include "main.h"
#include "utils.h"
#include "bwtindex.h"

/*"This BWT indexing program is part of the open-source BWA software that complies with the MIT license*/

#define PROGRAM_MACRO		"cushaw3 index"

bwt_t *bwt_pac2bwt(const char *fn_pac, int use_is);
void bwa_pac_rev_core(const char *fn, const char *fn_rev);

int cushaw2_index(int argc, char *argv[]) {
	int intVal;
	char *prefix = 0, *str, *str2, *str3;
	int c, is_color = 0;
	clock_t t;
	uint32_t refSeqLength = 0;
	bwt_t *bwt;
	int saInterval, powerOfInterval = 3;

	while ((c = getopt(argc, argv, "cp:i:")) >= 0) {
		switch (c) {
		case 'i':
			intVal = atoi(optarg);
			if (intVal < 0)
				intVal = 0;
			if (intVal > 6)
				intVal = 6;
			powerOfInterval = intVal;
			break;
		case 'p':
			prefix = strdup(optarg);
			break;
		case 'c':
			is_color = 1;
			break;
		default:
			return 1;
		}
	}
	if (optind + 1 > argc) {
		fprintf(stderr,
				"\nUsage: %s [options] <ref.fasta>\n", PROGRAM_MACRO);
		fprintf(stderr,
				"\t-p <string> (prefix of the index, default = fasta name)\n");
		fprintf(stderr,
				"\t-i <int> (interval for reduced suffix array [%d], 2^#INT)\n",
				powerOfInterval);
		fprintf(stderr, "\t-c (build color-space index)\n");
		return 1;
	}
	saInterval = pow(2, powerOfInterval);
	//fprintf(stderr, "SA interval: %d\n", saInterval);

	/*create file names*/
	if (prefix == 0){
		prefix = strdup(argv[optind]);
	}
	str = (char*) calloc(strlen(prefix) + 32, 1);
	str2 = (char*) calloc(strlen(prefix) + 32, 1);
	str3 = (char*) calloc(strlen(prefix) + 32, 1);

	/*pack the genome*/
	if (is_color == 0) { // nucleotide indexing
		gzFile fp = xzopen(argv[optind], "r");
		t = clock();
		fprintf(stderr, "[%s] Pack FASTA... ", PROGRAM_MACRO);
		bns_fasta2bntseq(fp, prefix);
		fprintf(stderr, "%.2f sec\n", (float) (clock() - t) / CLOCKS_PER_SEC);
		gzclose(fp);
	} else { // color indexing
		gzFile fp = xzopen(argv[optind], "r");
		strcat(strcpy(str, prefix), ".nt");
		t = clock();
		fprintf(stderr, "[%s] Pack nucleotide FASTA... ", PROGRAM_MACRO);
		bns_fasta2bntseq(fp, str);
		fprintf(stderr, "%.2f sec\n", (float) (clock() - t) / CLOCKS_PER_SEC);
		gzclose(fp);
		{
			char *tmp_argv[3];
			tmp_argv[0] = argv[0];
			tmp_argv[1] = str;
			tmp_argv[2] = prefix;
			t = clock();
			fprintf(stderr,
					"[%s] Convert nucleotide PAC to color PAC... ", PROGRAM_MACRO);
			bwa_pac2cspac(3, tmp_argv);
			fprintf(stderr, "%.2f sec\n",
					(float) (clock() - t) / CLOCKS_PER_SEC);
		}
	}

	/*reverse the packed genome*/
	strcpy(str, prefix);
	strcat(str, ".pac");
	strcpy(str2, prefix);
	strcat(str2, ".rpac");
	t = clock();
	fprintf(stderr, "[%s] Reverse the packed sequence... ", PROGRAM_MACRO);
	bwa_pac_rev_core(str, str2);
	fprintf(stderr, "%.2f sec\n", (float) (clock() - t) / CLOCKS_PER_SEC);

	/*reverse BWT*/
	strcpy(str, prefix);
	strcat(str, ".rpac");
	strcpy(str2, prefix);
	strcat(str2, ".rbwt");

	/*get the reference length*/
	refSeqLength = bwa_seq_len(str);
	fprintf(stderr, "[%s] Concatenated reference length: %u\n", PROGRAM_MACRO, refSeqLength);

	t = clock();
	fprintf(stderr,
			"[%s] Construct BWT for the reverse packed sequence...\n", PROGRAM_MACRO);
	if (refSeqLength > 4000000){/*more than 4 million*/
		/*using bwtsw type*/
		bwt_bwtgen(str, str2);
	}else{
		/*using is type*/
		bwt = bwt_pac2bwt(str, 1);
		bwt_dump_bwt(str2, bwt);
		bwt_destroy(bwt);
	}
	fprintf(stderr, "[%s] %.2f seconds elapse.\n", PROGRAM_MACRO,
			(float) (clock() - t) / CLOCKS_PER_SEC);

	/*update reverse BWT*/
	strcpy(str, prefix);
	strcat(str, ".rbwt");
	t = clock();
	fprintf(stderr, "[%s] Update reverse BWT... ", PROGRAM_MACRO);
	bwt = bwt_restore_bwt(str);
	bwt_bwtupdate_core(bwt);
	bwt_dump_bwt(str, bwt);
	bwt_destroy(bwt);
	fprintf(stderr, "%.2f sec\n", (float) (clock() - t) / CLOCKS_PER_SEC);


	/*reverse SA*/
	strcpy(str, prefix);
	strcat(str, ".rbwt");
	strcpy(str3, prefix);
	strcat(str3, ".rsa");
	t = clock();
	fprintf(stderr,
			"[%s] Construct SA from reverse BWT and Occ... ", PROGRAM_MACRO);
	bwt = bwt_restore_bwt(str);
	bwt_cal_sa(bwt, saInterval);
	bwt_dump_sa(str3, bwt);
	bwt_destroy(bwt);
	fprintf(stderr, "%.2f sec\n", (float) (clock() - t) / CLOCKS_PER_SEC);

	/*delete unused files*/
	strcpy(str, prefix);
	strcat(str, ".rpac");
	unlink(str);

	/*build a bitmap for unkown bases*/
	{
		gzFile fp = xzopen(argv[optind], "r");
		t = clock();
		fprintf(stderr, "[%s] Build a bitmap for unknown bases... ", PROGRAM_MACRO);
		bns_random_bitmap(fp, prefix, refSeqLength);
		fprintf(stderr, "%.2f sec\n", (float) (clock() - t) / CLOCKS_PER_SEC);
		gzclose(fp);
	}

	/*remove unused files*/
	free(str3);
	free(str2);
	free(str);
	free(prefix);
	return 0;
}
