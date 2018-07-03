#include "SingleEnd.h"
#include "PairedEnd.h"
#include "BSOptions.h"
#include "CSOptions.h"
extern "C" {
#include "bwtindex.h"
}

static int invokeAligner(Options* options, int argc, char* argv[]) {
	double stime, etime;

	/*get the startup time*/
	stime = Utils::getSysTime();

	/*parse the arguments*/
	if (!options->parse(argc, argv)) {
		options->printUsage();
		return 0;
	}

	/*create the genome object*/
	Genome* rgenome = new Genome(options);

	/*create the SAM output object*/
	SAM* sam = new SAM(options, rgenome);

	/*run the engines*/
	if (options->isPaired() == false) {
		/*for single-end alignment*/
		SingleEnd* single = new SingleEnd(options, rgenome, sam);

		/*run the alignment*/
		single->execute();

		/*release aligner*/
		delete single;
	} else {
		/*for paired-end alignment*/
		PairedEnd* paired = new PairedEnd(options, rgenome, sam);

		/*run the alignment*/
		paired->execute();

		/*release aligner*/
		delete paired;
	}

	/*release resources*/
	delete sam;
	delete rgenome;

	/*report the wall clock time*/
	etime = Utils::getSysTime();
	fprintf(stderr, "Overall time: %.2f seconds\n", etime - stime);

	return 0;
}
static void printCommands() {
	fprintf(stderr,
			"\n%s (v%s) is a fast gapped long-read aligner based on Burrows-Wheeler transform\n\n",
			PROGRAM_NAME, PROGRAM_VERSION);

	fprintf(stderr, "Usage:\t%s <command> [options]\n\n", PROGRAM_NAME);
	fprintf(stderr, "Command:\n");
	fprintf(stderr, "\tindex\t\tbuild BWT and FM-index\n");
	fprintf(stderr,
			"\talign\t\tperform base-space read alignments (e.g. Illumina, 454, Ion Torrent and PacBio)\n");
	fprintf(stderr,
			"\tcalign\t\tperform color-space read alignments (ABI SOLiD)\n");
	fprintf(stderr, "\n");
}
int main(int argc, char* argv[]) {
	Options* options = NULL;

	/*parse the command*/
	if (argc < 2) {
		printCommands();
		return 1;
	}

	if (!strcmp(argv[1], "index")) { /*invoke genome indexer*/
		cushaw2_index(argc - 1, argv + 1);
	} else if (!strcmp(argv[1], "align")) { /*base-space aligner*/
		options = new BSOptions();
		invokeAligner(options, argc - 1, argv + 1);
	} else if (!strcmp(argv[1], "calign")) {/*color-space aligner*/
		options = new CSOptions();
		invokeAligner(options, argc - 1, argv + 1);
	} else {
		Utils::log("Unknown command: %s\n", argv[1]);
		return 1;
	}
	/*release resources*/
	if (options) {
		delete options;
	}

	return 0;
}
