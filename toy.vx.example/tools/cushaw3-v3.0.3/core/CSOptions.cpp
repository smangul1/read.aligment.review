/*
 * CSOptions.cpp
 *
 *  Created on: Mar 19, 2013
 *      Author: yongchao
 */

#include "CSOptions.h"
/*
 * BSCSOptions.cpp
 *
 *  Created on: Mar 19, 2013
 *      Author: yongchao
 */

/*
 * CSOptions.cpp
 *
 *  Created on: Jan 10, 2012
 *      Author: yongchao
 */

#include "CSOptions.h"

CSOptions::CSOptions() :
		Options() {
	/*set default*/
	_setDefaults();
}

void CSOptions::printUsage() {
	fprintf(stderr, "\nUsage: %s calign -r bwt [options]\n", PROGRAM_NAME);
	fprintf(stderr, "\tbwt: the file name base for the reference genome\n");

	/*the file input options*/
	fprintf(stderr, "Input:\n");
	fprintf(stderr,
			"\t-r <string> (the file name base for the reference genome)\n");

	/*single-end reads*/
	fprintf(stderr,
			"\t-f <string> file1 [file2] (single-end sequence files in FASTA/FASTQ format)\n");
	fprintf(stderr,
			"\t-b <string> file1 [file2] (single-end sequence files in BAM format)\n");
	fprintf(stderr,
			"\t-s <string> file1 [file2] (single-end sequence files in SAM format)\n");

	/*paired-end/mate-paired reads*/
	fprintf(stderr,
			"\t-q <string> file1_1 file1_2  [file2_1 file2_2] (paired-end/mate-paired sequence files in FASTA/FASTQ format)\n");
	fprintf(stderr,
			"\t-mode <int> (paired-end/mate-paired reads, default = %d)\n",
			_pairingMode);
	fprintf(stderr, "\t\t0 means paired-end reads\n");
	fprintf(stderr, "\t\t1 means mate-paired reads\n");

	/*color-space trimming*/
	fprintf(stderr,
			"\t-trim_primer <int> (trim the primer base in color-space reads, default = %d)\n",
			_trimPrimer);

	/*output*/
	fprintf(stderr, "Output:\n");
	fprintf(stderr, "\t-o <string> (SAM output file path, default = %s)\n",
			_samFileName.length() > 0 ? _samFileName.c_str() : "STDOUT");
	fprintf(stderr,
			"\t-min_qual <int> (minimal mapping quality score of a reported alignment, default = %d)\n",
			_minMapQual);
	fprintf(stderr,
			"\t-multi <int> (output <= #int [1, %d] top alignments, default = %d)\n",
			MAX_MULTI_ALIGN, _maxMultiAligns);

	/*read group information*/
	fprintf(stderr,
			"\n\t/*read group information in SAM header as follows*/\n");
	fprintf(stderr, "\t-rgid <string> (read group identifier [tag RG:ID])\n");
	fprintf(stderr,
			"\t-rgsm <string> (read group sample name [tag RG:SM], required if #rgid is given)\n");
	fprintf(stderr,
			"\t-rglb <string> (read group library [tag RG:LD], ineffective if #rgid is not given)\n");
	fprintf(stderr,
			"\t-rgpl <string> (read group platform/technology [tag RG:PL], ineffective if #rigid is not given)\n\t\tsupported values: capillary, ls454, illumina, solid, helicos, iontorrent, and pacbio\n");
	fprintf(stderr,
			"\t-rgpu <string> (read group platform unit identifier [tag RG:PU], ineffective if #rgid is not given)\n");
	fprintf(stderr,
			"\t-rgcn <string> (name of sequencing center produced the read [tag RG:CN], ineffiective if #rgid is not given)\n");
	fprintf(stderr,
			"\t-rgds <string> (description about the reads [tag RG:DS], ineffiective if #rgid is not given)\n");
	fprintf(stderr,
			"\t-rgdt <string> (date on which the run was produced [tag RG:DT], ineffiective if #rgid is not given)\n");
	fprintf(stderr,
			"\t-rgpi <string> (predicated median insert size [tag RG:PI], ineffiective if #rgid is not given)\n");

	/*scoring scheme*/
	fprintf(stderr, "Scoring:\n");
	fprintf(stderr, "\t-match <int> (score for a match, default = %d)\n",
			_match);
	fprintf(stderr,
			"\t-mismatch <int> (penalty for a mismatch, default = %d)\n",
			_mismatch);
	fprintf(stderr, "\t-gopen <int> (gap open penalty, default = %d)\n",
			_gapOpen);
	fprintf(stderr, "\t-gext <int> (gap extension penalty, default = %d)\n",
			_gapExtend);

	/*alignment options*/
	fprintf(stderr, "Align:\n");
	fprintf(stderr, "\t-atype <int> (alignment type used, default = %d)\n",
			_alignType);
	fprintf(stderr, "\t\t0 [local alignment]\n");
	fprintf(stderr, "\t\t1 [local + semiglobal alignment]\n");

	fprintf(stderr,
			"\t-mask_amb <int> (masking ambiguous bases in the reference, default = %d)\n",
			_maskAmbiguous);
	fprintf(stderr,
			"\t-min_score <int> (minimal optimal local alignment score divided by matching score, default = %d)\n",
			_minAlginScore);
	fprintf(stderr,
			"\t-min_id <float> (minimal identity of optical local alignments, default = %.2f)\n",
			_minIdentity);
	fprintf(stderr,
			"\t-gmin_id <float> (minimal identity of optimal semiglobal alignments, default = %.2f)\n",
			_gminIdentity);

	fprintf(stderr,
			"\t-min_ratio <float> (minimal ratio of reads in optimal local/semiglobal alignments, default = %.2f)\n",
			_minRatio);
	fprintf(stderr,
			"\t-max_edit_dist <int> (maximum edit distance in optimal local/semiglobal alignments, default = %d [-1 means disabled])\n",
			_maxEditDistance);

	/*seeding options*/
	fprintf(stderr, "Seed:\n");
	fprintf(stderr,
			"\t-min_seed <int> (lower bound of minimal seed size, default = %d)\n",
			_lowerMinSeedSize);
	fprintf(stderr,
			"\t-max_seed <int> (upper bound of minimal seed size, default = %d)\n",
			_upperMinSeedSize);
	fprintf(stderr,
			"\t-miss_prob <float> (missing probability to estimate the seed sizes, default = %.2f)\n",
			_missProb);
	fprintf(stderr,
			"\t-max_occ <int> (maximal number of occurrences per seed, default = %d)\n",
			_maxSeedOcc);
	fprintf(stderr,
			"\t-max_seeds_per_batch <int> (maximum number of seeds per batch for top hit selection, default = %d)\n",
			_maxSeedsPerBatch);

	/*Mate-paired alignment*/
	fprintf(stderr, "Pairing:\n");
	fprintf(stderr,
			"\t-sensitive_pairing <int> (more sensitive read pairing, default = %d)\n",
			_sensitivePairing);
	if (_estInsertSize) {
		fprintf(stderr,
				"\t-avg_ins <int> (insert size for paired-end/mate-paired reads, estimated from input if not specified)\n");

		fprintf(stderr,
				"\t-ins_std <int> (standard deviation of insert size, estimated from input if not specified)\n");
		fprintf(stderr,
				"\t-ins_npairs <int> (top number of read pairs for insert size estimation [#int times 0x%x], default = %d)\n",
				INS_SIZE_EST_MULTIPLE, _topReadsEstIns / INS_SIZE_EST_MULTIPLE);
		fprintf(stderr,
				"\t-ins_mapq <int> (minimal mapping quality score of a SE alignment for insert size estimation, default = %d)\n",
				_mapQualReliable);
	} else {
		fprintf(stderr,
				"\t-avg_ins <int> (insert size for paired-end/mate-paired reads, default = %d)\n",
				_insertSize);

		fprintf(stderr,
				"\t-ins_std <int> (standard deviation of insert size for paired-end/mate-paired reads, default = %d)\n",
				_stdInsertsize);
	}
	fprintf(stderr,
			"\t-max_seedpairs <int> (maximal number of seed pairs per read pair, default = %d)\n",
			_maxSeedPairs);
	fprintf(stderr,
			"\t-no_rescue (do not rescue the mate using Smith-Waterman for an un-paired read)\n");
	fprintf(stderr,
			"\t-rescue_top_hits <int> (attempt the second rescuing from #int top hits for each end, default = %d)\n",
			_rescueTwice);

	/*compute*/
	fprintf(stderr, "Compute:\n");
	fprintf(stderr, "\t-t <int> (number of threads, default = %d)\n",
			_numThreads);
	fprintf(stderr, "Others:\n");
	fprintf(stderr, "\t-h <print out the usage of the program)\n");
}
void CSOptions::_setDefaults() {
	/*parameters for file input*/
	_bwtFileBase = "";

	/*empty string means outputing to STDOUT*/
#ifdef COMPRESSED_SAM
	_samFileName="/dev/stdout";
#else
	_samFileName = "";
#endif
	_readsFileNames.clear();

	/*parameters for alignment*/
	_minRatio = DEFAULT_SEQ_LENGTH_RATIO; /*minimal portion for bases in a short read*/
	_minIdentity = DEFAULT_MIN_IDENTITY; /*minimal identity*/
	_gminIdentity = DEFAULT_GLOBAL_MIN_IDENTITY; /*minimal identiy for semiglobal alignment*/
	_numThreads = 1; /*the number of threads used for the alignment*/
	_estInsertSize = true; /*estimate the insert size from input paired-rend reads*/
	_insertSize = 500; /*the insert size for paired-end/mate-paired reads*/
	_stdInsertsize = (int) (_insertSize * 0.1); /*the standard deviation of insert size for paired-end/mate-paired reads*/
	_topReadsEstIns = 0x10000;
	_mapQualReliable = 20; /*for insert-size estimation*/
	_minMapQual = 0; /*for alignment output*/
	_maxMultiAligns = 10;
	_missProb = 0.04;
	_paired = false;
	_sensitivePairing = false; /*for paired-end alignment*/
	_maskAmbiguous = true;
	_pairingMode = MATE_PAIRED_READS;
	_lowerMinSeedSize = 13;
	_upperMinSeedSize = GLOBAL_MAX_SEED_SIZE;
	_maxSeedOcc = GLOBAL_MAX_NUM_SEED_REPEATS;
	_maxSeedPairs = GLOBAL_MAX_NUM_SEED_PAIRS;
	_rescueMate = true;
	_rescueTwice = 100;
	_maxSeedsPerBatch = 16384;

	/*color-space*/
	_colorspace = true;
	_trimPrimer = true;
	_maxEditDistance = -1; /*maximum edit distance in the alignment for color-space reads*/

	/*penalties*/
	_match = 1; /*score for matching*/
	_mismatch = 3; /*penalty for mismatching*/
	_gapOpen = 5; /*penalty for gap open*/
	_gapExtend = 2; /*penalty for gap extension*/
	_minAlginScore = 10 * _match; /*minial local alignment score*/
	_alignType = HAVE_LOCAL_SEMI_GLOBAL_ALIGN; /*using local + semiglobal alignment*/

	/*estimate errors*/
	_numErrorTable = new int[MAX_SEQ_LENGTH + 1];
	if (_numErrorTable == NULL) {
		Utils::exit("Memory allocation failed in function %s line %d\n",
				__FUNCTION__, __LINE__);
	}
	_minSeedSizeTable = new int[MAX_SEQ_LENGTH + 1];
	if (_minSeedSizeTable == NULL) {
		Utils::exit("Memory allocation failed in function %s line %d\n",
				__FUNCTION__, __LINE__);
	}

	/*global lock*/
	pthread_mutex_init(&globalMutex, NULL);
}
bool CSOptions::parse(int argc, char* argv[]) {
	bool newInsertSize = false;
	bool newStdInsertSize = false;
	int intVal;
	double floatVal;
	int argind = 1;

	if (argc < 2) {
		return false;
	}

	/*check the help*/
	if (!strcmp(argv[argind], "-h") || !strcmp(argv[argind], "-?")) {
		return false;
	}

	/*print out the command line*/
	fprintf(stderr, "Command: ");
	for (int i = 0; i < argc; ++i) {
		fprintf(stderr, "%s ", argv[i]);
	}
	fputc('\n', stderr);

	/*for the other options*/
	bool first = true;
	bool done = false;
	while (argind < argc) {
		/*single-end sequence files*/
		if (!strcmp(argv[argind], "-r")) {
			argind++;
			if (argind < argc) {
				_bwtFileBase = argv[argind];
				argind++;
			} else {
				Utils::log("not specify value for the parameter %s\n",
						argv[argind - 1]);
				return false;
			}

		} else if (!strcmp(argv[argind], "-f")) {
			/*specify the inputs are single-end*/
			if (first) {
				first = false;
				setPaired(false);
			} else if (isPaired()) {
				Utils::log(
						"Cannot specify single-end and paired-end/mate-paired reads together\n");
				return false;
			}

			/*increase the argument index*/
			argind++;
			done = false;
			while (argind < argc
					&& (argv[argind][0] != '-'
							|| (argv[argind][0] == '-'
									&& argv[argind][1] == '\0'))) {
				//test file
				if (strcmp(argv[argind], "-") && !Utils::exists(argv[argind])) {
					Utils::log("The file %s does not exist\n", argv[argind]);
					return false;
				}
				//get the file name
				_readsFileNames.push_back(
						make_pair(string(argv[argind]), FILE_TYPE_FASTX));

				//increase the argument index
				argind++;
				done = true;
			}
			if (done == false) {
				Utils::log("At least one input file must be specified\n");
				return false;
			}
		} else if (!strcmp(argv[argind], "-b")) {
			/*specify the inputs are single-end*/
			if (first) {
				first = false;
				setPaired(false);
			} else if (isPaired()) {
				Utils::log(
						"Cannot specify single-end and paired-end/mate-paired reads together\n");
				return false;
			}

			/*increase the argument index*/
			argind++;
			done = false;
			while (argind < argc
					&& (argv[argind][0] != '-'
							|| (argv[argind][0] == '-'
									&& argv[argind][1] == '\0'))) {
				//test file
				if (strcmp(argv[argind], "-") && !Utils::exists(argv[argind])) {
					Utils::log("The file %s does not exist\n", argv[argind]);
					return false;
				}
				//get the file name
				_readsFileNames.push_back(
						make_pair(string(argv[argind]), FILE_TYPE_BAM));

				//increase the argument index
				argind++;
				done = true;
			}
			if (done == false) {
				Utils::log("At least one input file must be specified\n");
				return false;
			}
		} else if (!strcmp(argv[argind], "-s")) {
			/*specify the inputs are single-end*/
			if (first) {
				first = false;
				setPaired(false);
			} else if (isPaired()) {
				Utils::log(
						"Cannot specify single-end and paired-end/mate-paired reads together\n");
				return false;
			}

			/*increase the argument index*/
			argind++;
			done = false;
			while (argind < argc
					&& (argv[argind][0] != '-'
							|| (argv[argind][0] == '-'
									&& argv[argind][1] == '\0'))) {
				//test file
				if (strcmp(argv[argind], "-") && !Utils::exists(argv[argind])) {
					Utils::log("The file %s does not exist\n", argv[argind]);
					return false;
				}
				//get the file name
				_readsFileNames.push_back(
						make_pair(string(argv[argind]), FILE_TYPE_SAM));

				//increase the argument index
				argind++;
				done = true;
			}
			if (done == false) {
				Utils::log("At least one input file must be specified\n");
				return false;
			}

			//input files are paired in FASTA format
		} else if (!strcmp(argv[argind], "-q")) {
			if (first) {
				first = false;
				setPaired(1);
			} else if (!isPaired()) {
				Utils::log(
						"Cannot specify single-end and paired-end/mate-paired reads together\n");
				return false;
			}

			++argind;
			done = false;
			if (argind + 1 < argc && argv[argind][0] != '-'
					&& (argv[argind + 1][0] != '-'
							|| (argv[argind + 1][0] == '-'
									&& argv[argind + 1][1] == '\0'))) {
				//save the two files
				//test file
				if (strcmp(argv[argind], "-") && !Utils::exists(argv[argind])) {
					Utils::log("The file %s does not exist\n", argv[argind]);
					return false;
				}
				_readsFileNames.push_back(
						make_pair(string(argv[argind]), FILE_TYPE_FASTX));
				argind++;

				//test file
				if (strcmp(argv[argind], "-") && !Utils::exists(argv[argind])) {
					Utils::log("The file %s does not exist\n", argv[argind]);
					return false;
				}
				_readsFileNames.push_back(
						make_pair(string(argv[argind]), FILE_TYPE_FASTX));
				argind++;

				done = true;
			}
			if (!done) {
				Utils::log(
						"Two paired input files should be specified for -fq\n");
				return false;
			}
		} else if (!strcmp(argv[argind], "-trim_primer")) {
			intVal = 1;
			argind++;
			if (argind < argc) {
				sscanf(argv[argind], "%d", &intVal);
				argind++;
			} else {
				Utils::log("not specify value for the parameter %s\n",
						argv[argind - 1]);
				return false;
			}
			_trimPrimer = intVal != 0 ? true : false;
		} else if (!strcmp(argv[argind], "-mode")) {
			intVal = 1;
			argind++;
			if (argind < argc) {
				sscanf(argv[argind], "%d", &intVal);
				argind++;
			} else {
				Utils::log("not specify value for the parameter %s\n",
						argv[argind - 1]);
				return false;
			}
			_pairingMode = (intVal == 0) ? PAIRED_END_READS : MATE_PAIRED_READS;
		} else if (!strcmp(argv[argind], "-o")) {
			argind++;
			if (argind < argc) {
				if (strlen(argv[argind]) > 0) {
					_samFileName = argv[argind];
					argind++;
				}
			} else {
				Utils::log("not specify value for the parameter %s\n",
						argv[argind - 1]);
				return false;
			}
		} else if (!strcmp(argv[argind], "-multi")) {
			intVal = 1;
			argind++;
			if (argind < argc) {
				sscanf(argv[argind], "%d", &intVal);
				if (intVal < 1) {
					intVal = 1;
				}
				if (intVal > MAX_MULTI_ALIGN) {
					intVal = MAX_MULTI_ALIGN;
				}
				_maxMultiAligns = intVal;
				argind++;
			} else {
				Utils::log("not specify value for the parameter %s\n",
						argv[argind - 1]);
				return false;
			}
		} else if (!strcmp(argv[argind], "-min_qual")) {
			intVal = 1;
			argind++;
			if (argind < argc) {
				sscanf(argv[argind], "%d", &intVal);
				if (intVal < 0) {
					intVal = 0;
				}
				_minMapQual = intVal;
				argind++;
			} else {
				Utils::log("not specify value for the parameter %s\n",
						argv[argind - 1]);
				return false;
			}
			/*read group information*/
		} else if (!strcmp(argv[argind], "-rgid")) { /*read group identifier*/
			argind++;
			if (argind < argc) {
				if (strlen(argv[argind]) > 0) {
					_rgID = argv[argind];
					argind++;
				}
			} else {
				Utils::log("not specify value for the parameter %s\n",
						argv[argind - 1]);
				return false;
			}
		} else if (!strcmp(argv[argind], "-rgsm")) { /*read group sample name*/
			argind++;
			if (argind < argc) {
				if (strlen(argv[argind]) > 0) {
					_rgSM = argv[argind];
					argind++;
				}
			} else {
				Utils::log("not specify value for the parameter %s\n",
						argv[argind - 1]);
				return false;
			}
		} else if (!strcmp(argv[argind], "-rglb")) { /*read group library*/
			argind++;
			if (argind < argc) {
				if (strlen(argv[argind]) > 0) {
					_rgLB = argv[argind];
					argind++;
				}
			} else {
				Utils::log("not specify value for the parameter %s\n",
						argv[argind - 1]);
				return false;
			}
		} else if (!strcmp(argv[argind], "-rgpl")) { /*read group platform*/
			argind++;
			if (argind < argc) {
				if (strlen(argv[argind]) > 0) {
					_rgPL = argv[argind];
					argind++;
				}
			} else {
				Utils::log("not specify value for the parameter %s\n",
						argv[argind - 1]);
				return false;
			}
		} else if (!strcmp(argv[argind], "-rgpu")) { /*read group platform unit*/
			argind++;
			if (argind < argc) {
				if (strlen(argv[argind]) > 0) {
					_rgPU = argv[argind];
					argind++;
				}
			} else {
				Utils::log("not specify value for the parameter %s\n",
						argv[argind - 1]);
				return false;
			}
		} else if (!strcmp(argv[argind], "-rgcn")) { /*name of sequencing center producing the read*/
			argind++;
			if (argind < argc) {
				if (strlen(argv[argind]) > 0) {
					_rgCN = argv[argind];
					argind++;
				}
			} else {
				Utils::log("not specify value for the parameter %s\n",
						argv[argind - 1]);
				return false;
			}
		} else if (!strcmp(argv[argind], "-rgds")) { /*description about the reads*/
			argind++;
			if (argind < argc) {
				if (strlen(argv[argind]) > 0) {
					_rgDS = argv[argind];
					argind++;
				}
			} else {
				Utils::log("not specify value for the parameter %s\n",
						argv[argind - 1]);
				return false;
			}
		} else if (!strcmp(argv[argind], "-rgdt")) { /*date the run was produced*/
			argind++;
			if (argind < argc) {
				if (strlen(argv[argind]) > 0) {
					_rgDT = argv[argind];
					argind++;
				}
			} else {
				Utils::log("not specify value for the parameter %s\n",
						argv[argind - 1]);
				return false;
			}
		} else if (!strcmp(argv[argind], "-rgpi")) { /*predicated median insert size*/
			argind++;
			if (argind < argc) {
				if (strlen(argv[argind]) > 0) {
					_rgPI = argv[argind];
					argind++;
				}
			} else {
				Utils::log("not specify value for the parameter %s\n",
						argv[argind - 1]);
				return false;
			}
		} else if (!strcmp(argv[argind], "-min_ratio")) {
			floatVal = 1;
			argind++;
			if (argind < argc) {
				sscanf(argv[argind], "%lf", &floatVal);
				if (floatVal < 0) {
					floatVal = 0;
				}
				if (floatVal > 1) {
					floatVal = 1;
				}
				argind++;
				_minRatio = floatVal;
			} else {
				Utils::log("not specify value for the parameter %s\n",
						argv[argind - 1]);
				return false;
			}
		} else if (!strcmp(argv[argind], "-min_id")) { /*minimal identity*/
			floatVal = 1;
			argind++;
			if (argind < argc) {
				sscanf(argv[argind], "%lf", &floatVal);
				if (floatVal < 0) {
					floatVal = 0;
				}
				if (floatVal > 1) {
					floatVal = 1;
				}
				argind++;
			} else {
				Utils::log("not specify value for the parameter %s\n",
						argv[argind - 1]);
				return false;
			}
			_minIdentity = floatVal;
		} else if (!strcmp(argv[argind], "-gmin_id")) { /*minimal identity*/
			floatVal = 1;
			argind++;
			if (argind < argc) {
				sscanf(argv[argind], "%lf", &floatVal);
				if (floatVal < 0) {
					floatVal = 0;
				}
				if (floatVal > 1) {
					floatVal = 1;
				}
				argind++;
			} else {
				Utils::log("not specify value for the parameter %s\n",
						argv[argind - 1]);
				return false;
			}
			_gminIdentity = floatVal;
		} else if (!strcmp(argv[argind], "-max_edit_dist")) { /*minimal identity*/
			intVal = 1;
			argind++;
			if (argind < argc) {
				sscanf(argv[argind], "%d", &intVal);
				if (intVal < 0) {
					intVal = -1;
				}
				_maxEditDistance = intVal;
				argind++;
			} else {
				Utils::log("not specify value for the parameter %s\n",
						argv[argind - 1]);
				return false;
			}
		} else if (!strcmp(argv[argind], "-min_score")) {
			intVal = 1;
			argind++;
			if (argind < argc) {
				sscanf(argv[argind], "%d", &intVal);
				if (intVal < 1) {
					intVal = 1;
				}
				_minAlginScore = intVal;
				argind++;
			} else {
				Utils::log("not specify value for the parameter %s\n",
						argv[argind - 1]);
				return false;
			}
		} else if (!strcmp(argv[argind], "-mask_amb")) {
			intVal = 1;
			argind++;
			if (argind < argc) {
				sscanf(argv[argind], "%d", &intVal);
				_maskAmbiguous = intVal ? true : false;
				argind++;
			} else {
				Utils::log("not specify value for the parameter %s\n",
						argv[argind - 1]);
				return false;
			}
		} else if (!strcmp(argv[argind], "-atype")) {
			intVal = HAVE_LOCAL_ALIGN;
			argind++;
			if (argind < argc) {
				sscanf(argv[argind], "%d", &intVal);
				if (intVal >= 0
						&& intVal < sizeof(_atypes) / sizeof(_atypes[0])) {
					_alignType = _atypes[intVal];
				} else {
					Utils::log("The specified value (%d) is not supported\n",
							intVal);
					return false;
				}
				argind++;
			} else {
				Utils::log("not specify value for the parameter %s\n",
						argv[argind - 1]);
				return false;
			}
		} else if (!strcmp(argv[argind], "-miss_prob")) {
			floatVal = 1;
			argind++;
			if (argind < argc) {
				sscanf(argv[argind], "%lf", &floatVal);
				if (floatVal < 0)
					floatVal = 0;
				if (floatVal > 1)
					floatVal = 1;
				argind++;
				_missProb = floatVal;
			} else {
				Utils::log("not specify value for the parameter %s\n",
						argv[argind - 1]);
				return false;
			}

		} else if (!strcmp(argv[argind], "-min_seed")) {
			intVal = 1;
			argind++;
			if (argind < argc) {
				sscanf(argv[argind], "%d", &intVal);
				if (intVal < GLOBAL_MIN_SEED_SIZE) {
					intVal = GLOBAL_MIN_SEED_SIZE;
				} else if (intVal > GLOBAL_MAX_SEED_SIZE) {
					intVal = GLOBAL_MAX_SEED_SIZE;
				}
				_lowerMinSeedSize = intVal;
				argind++;
			} else {
				Utils::log("not specify value for the parameter %s\n",
						argv[argind - 1]);
				return false;
			}
		} else if (!strcmp(argv[argind], "-max_seed")) {
			intVal = 1;
			argind++;
			if (argind < argc) {
				sscanf(argv[argind], "%d", &intVal);
				if (intVal < GLOBAL_MIN_SEED_SIZE) {
					intVal = GLOBAL_MIN_SEED_SIZE;
				} else if (intVal > GLOBAL_MAX_SEED_SIZE) {
					intVal = GLOBAL_MAX_SEED_SIZE;
				}
				_upperMinSeedSize = intVal;
				argind++;
			} else {
				Utils::log("not specify value for the parameter %s\n",
						argv[argind - 1]);
				return false;
			}
		} else if (!strcmp(argv[argind], "-max_occ")) {
			intVal = 1;
			argind++;
			if (argind < argc) {
				sscanf(argv[argind], "%d", &intVal);
				if (intVal < 1) {
					intVal = 1;
				}
				_maxSeedOcc = intVal;
				argind++;
			} else {
				Utils::log("not specify value for the parameter %s\n",
						argv[argind - 1]);
				return false;
			}
		} else if (!strcmp(argv[argind], "-max_seeds_per_batch")) {
			intVal = 1;
			argind++;
			if (argind < argc) {
				sscanf(argv[argind], "%d", &intVal);
				if (intVal < 1024) {
					intVal = 1024;
				}
				_maxSeedsPerBatch = intVal;
				argind++;
			} else {
				Utils::log("not specify value for the parameter %s\n",
						argv[argind - 1]);
				return false;
			}
		} else if (!strcmp(argv[argind], "-max_seed")) {
			intVal = 1;
			argind++;
			if (argind < argc) {
				sscanf(argv[argind], "%d", &intVal);
				if (intVal < GLOBAL_MIN_SEED_SIZE) {
					intVal = GLOBAL_MIN_SEED_SIZE;
				} else if (intVal > GLOBAL_MAX_SEED_SIZE) {
					intVal = GLOBAL_MAX_SEED_SIZE;
				}
				_upperMinSeedSize = intVal;
				argind++;
			} else {
				Utils::log("not specify value for the parameter %s\n",
						argv[argind - 1]);
				return false;
			}
		} else if (!strcmp(argv[argind], "-match")) {
			intVal = 1;
			argind++;
			if (argind < argc) {
				sscanf(argv[argind], "%d", &intVal);
				if (intVal < 1)
					intVal = 1;

				argind++;
				_match = intVal;
			} else {
				Utils::log("not specify value for the parameter %s\n",
						argv[argind - 1]);
				return false;
			}
		} else if (!strcmp(argv[argind], "-mismatch")) {
			intVal = 1;
			argind++;
			if (argind < argc) {
				sscanf(argv[argind], "%d", &intVal);
				if (intVal < 0)
					intVal = 0;

				argind++;
				_mismatch = intVal;
			} else {
				Utils::log("not specify value for the parameter %s\n",
						argv[argind - 1]);
				return false;
			}
		} else if (!strcmp(argv[argind], "-gopen")) {
			intVal = 1;
			argind++;
			if (argind < argc) {
				sscanf(argv[argind], "%d", &intVal);
				if (intVal < 0)
					intVal = 0;

				argind++;
				_gapOpen = intVal;
			} else {
				Utils::log("not specify value for the parameter %s\n",
						argv[argind - 1]);
				return false;
			}
		} else if (!strcmp(argv[argind], "-gext")) {
			intVal = 1;
			argind++;
			if (argind < argc) {
				sscanf(argv[argind], "%d", &intVal);
				if (intVal < 0)
					intVal = 0;

				argind++;
				_gapExtend = intVal;
			} else {
				Utils::log("not specify value for the parameter %s\n",
						argv[argind - 1]);
				return false;
			}
		} else if (!strcmp(argv[argind], "-sensitive_pair")) {
			intVal = 1;
			argind++;
			if (argind < argc) {
				sscanf(argv[argind], "%d", &intVal);

				argind++;
				_sensitivePairing = intVal ? true : false;
			} else {
				Utils::log("not specify value for the parameter %s\n",
						argv[argind - 1]);
				return false;
			}
		} else if (!strcmp(argv[argind], "-avg_ins")) {
			intVal = 1;
			argind++;
			if (argind < argc) {
				sscanf(argv[argind], "%d", &intVal);
				if (intVal < 0)
					intVal = 0;

				argind++;
				_insertSize = intVal;
				newInsertSize = true;
			} else {
				Utils::log("not specify value for the parameter %s\n",
						argv[argind - 1]);
				return false;
			}
		} else if (!strcmp(argv[argind], "-ins_std")) {
			intVal = 1;
			argind++;
			if (argind < argc) {
				sscanf(argv[argind], "%d", &intVal);
				if (intVal < 0)
					intVal = 0;

				argind++;
				_stdInsertsize = intVal;
				newStdInsertSize = true;
			} else {
				Utils::log("not specify value for the parameter %s\n",
						argv[argind - 1]);
				return false;
			}
		} else if (!strcmp(argv[argind], "-ins_npairs")) {
			intVal = 1;
			argind++;
			if (argind < argc) {
				sscanf(argv[argind], "%d", &intVal);
				if (intVal < 1)
					intVal = 1;

				argind++;
				_topReadsEstIns = intVal * INS_SIZE_EST_MULTIPLE;
			} else {
				Utils::log("not specify value for the parameter %s\n",
						argv[argind - 1]);
				return false;
			}
		} else if (!strcmp(argv[argind], "-ins_mapq")) {
			intVal = 1;
			argind++;
			if (argind < argc) {
				sscanf(argv[argind], "%d", &intVal);
				if (intVal < 0)
					intVal = 0;

				if (intVal > 255) {
					intVal = 255;
				}

				argind++;
				_mapQualReliable = intVal;
			} else {
				Utils::log("not specify value for the parameter %s\n",
						argv[argind - 1]);
				return false;
			}

		} else if (!strcmp(argv[argind], "-no_rescue")) {
			++argind;
			_rescueMate = false;
		} else if (!strcmp(argv[argind], "-rescue_top_hits")) {
			intVal = 1;
			argind++;
			if (argind < argc) {
				sscanf(argv[argind], "%d", &intVal);
				if (intVal < 0)
					intVal = 0;

				argind++;
				_rescueTwice = intVal;
			} else {
				Utils::log("not specify value for the parameter %s\n",
						argv[argind - 1]);
				return false;
			}
		} else if (!strcmp(argv[argind], "-t")) {
			intVal = 1;
			argind++;
			if (argind < argc) {
				sscanf(argv[argind], "%d", &intVal);
				if (intVal < 1)
					intVal = 1;

				argind++;
				_numThreads = intVal;
			} else {
				Utils::log("not specify value for the parameter %s\n",
						argv[argind - 1]);
				return false;
			}
		} else if (!strcmp(argv[argind], "-h")) {
			return false;
		} else {
			Utils::log("Unknown parameter: %s\n", argv[argind]);
			return false;
		}
	}
	/*check the genome*/
	if (_bwtFileBase.size() == 0) {
		Utils::log(
				"The reference genome must be specified using option \"-r\"");
		return false;
	}

	/*re-check the insert size*/
	if (newInsertSize && !newStdInsertSize) {
		/*if the standard deviation is not specified, just use 10% of the insert size*/
		_stdInsertsize = (int) (_insertSize * 0.1);
	}

	/*check if insert size is etimated automatically*/
	if (!newInsertSize && !newStdInsertSize) {
		_estInsertSize = true;
	} else {
		_estInsertSize = false;
	}

	/*check the seed size*/
	if (_upperMinSeedSize < _lowerMinSeedSize) {
		Utils::exit(
				"The upper bound of the minimal seed size (%d) is not allowed to be less than its lower bound (%d)\n",
				_upperMinSeedSize, _lowerMinSeedSize);
	}
	/*check the read group header information*/
	if (_rgID.size() > 0) {
		if (_rgSM.size() == 0) {
			Utils::exit(
					"Must specify read group sample name if read group identifier is specified\n");
		}
	}

	/*update the minimal local alignment score*/
	_minAlginScore *= _match;

	/*initialize the table*/
	int seedSize, numErrors;
	for (uint32_t length = 1; length <= MAX_SEQ_LENGTH; ++length) {
		numErrors = estimateNumErrors(length, _missProb);
		_numErrorTable[length] = numErrors;

		/*estimate the minimal seed size according to dove hole principle*/
		seedSize = length / (numErrors + 1);
		seedSize = max(seedSize, _lowerMinSeedSize);
		seedSize = min(seedSize, _upperMinSeedSize);

		/*save the minimal seed Size*/
		_minSeedSizeTable[length] = seedSize;
	}

	/*check the maximum edit distance*/
	if (_maxEditDistance < 0) {
		_maxEditDistance = INT_MAX;
	}

	return true;
}

