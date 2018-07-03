/*
 * SAM.cpp
 *
 *  Created on: Jan 5, 2012
 *      Author: yongchao
 */
#include "SAM.h"
#include "Utils.h"

SAM::SAM(Options* options, Genome* genome) {
	_fileName = options->getSamFileName();
	_numThreads = options->getNumThreads();
	_file = NULL;
	_genome = genome;
	_colorspace = options->isColorSpace();
	_pairingMode = options->getPairingMode();

	/*open the file*/
	open();

	/*print out the header*/
	fprintf(_file, "@HD\tVN:1.0\tSO:unsorted\n");

	/*print out the genome sequence information in SAM format*/
	_genome->genomeNamesOut(_file);

	/*print out the read group header information*/
	_groupID = NULL;
	_groupLB = NULL;
	if (options->_rgID.length() > 0) {
		_groupID = (char*) options->_rgID.c_str();
		/*tag ID*/
		fprintf(_file, "@RG\tID:%s\tSM:%s", options->_rgID.c_str(),
				options->_rgSM.c_str());
		/*tag LB*/
		if (options->_rgLB.length() > 0) {
			fprintf(_file, "\tLB:%s", options->_rgLB.c_str());
			_groupLB = (char*)options->_rgLB.c_str();
		}
		/*tag PL*/
		if (options->_rgPL.length() > 0) {
			fprintf(_file, "\tPL:%s", options->_rgPL.c_str());
		}
		/*tag PU*/
		if (options->_rgPU.length() > 0) {
			fprintf(_file, "\tPU:%s", options->_rgPU.c_str());
		}
		/*tag CN*/
		if (options->_rgCN.length() > 0) {
			fprintf(_file, "\tCN:%s", options->_rgCN.c_str());
		}
		/*tag DS*/
		if (options->_rgDS.length() > 0) {
			fprintf(_file, "\tDS:%s", options->_rgDS.c_str());
		}
		/*tag DT*/
		if (options->_rgDT.length() > 0) {
			fprintf(_file, "\tDT:%s", options->_rgDT.c_str());
		}
		/*tag PI*/
		if (options->_rgPI.length() > 0) {
			fprintf(_file, "\tPI:%s", options->_rgPI.c_str());
		}
		/*end of the line*/
		fputc('\n', _file);
	}

	/*print out the aligner information*/
	fprintf(_file, "@PG\tID:%s\tVN:%s\n", PROGRAM_NAME, PROGRAM_VERSION);
}
SAM::~SAM() {
	/*close the file*/
	close();

}
void SAM::open() {
	/*if the file is open, close it*/

	if (_file) {
		close();
	}

	/*re-open it*/
	if (_fileName.length() > 0) {
		/*open the output file*/
		_file = fopen(_fileName.c_str(), "wb");
	} else {
		/*open the standard output*/
		_file = stdout;
	}

	if (_file == NULL) {
		Utils::exit("Failed to open file %s in function %s line %d\n",
				_fileName.c_str(), __FUNCTION__, __LINE__);
	}
}
void SAM::close() {

	if (_file) {
		fclose(_file);
	}
	_file = NULL;

}
/*print the unaligned read information*/
void SAM::print(Sequence& seq, int _flags) {
	int flags = SAM_FSU | _flags;

	/*print out query name, bitwise-flag, and reference sequence name*/
	fprintf(_file, "%s\t%d\t*\t", seq._name, flags);

	//print 1-based leftmost mapping position, mapping quality (phred-scale)*/
	fprintf(_file, "0\t0\t");

	//print extended CIGAR
	fputs("*\t", _file);

	//print paired-end information, INAVAILABLE
	fprintf(_file, "*\t0\t0\t");

	//print the query sequence
	uint8_t* bases = seq._bases;
	for (uint32_t i = 0; i < seq._length; ++i) {
		fputc(decode(bases[i]), _file);
	}
	fputc('\t', _file);

	//print the quality scores if available
	if (seq._quals) {
		uint8_t* quals = seq._quals;
		for (uint32_t i = 0; i < seq._length; ++i) {
			fputc(*quals, _file);
			++quals;
		}
	} else {
		fputc('*', _file);
	}
	/*print tags*/
	fprintf(_file, "\tPG:Z:%s", PROGRAM_NAME);
	if (_groupID) {
		fprintf(_file, "\tRG:Z:%s", _groupID);
	}
	if(_groupLB){
		fprintf(_file, "\tLB:Z:%s", _groupLB);
	}
	fputc('\n', _file);
}
/*print the alignment information through the dynamic programming*/
void SAM::print(Sequence& seq, Mapping& mapping, int _flags) {
	int32_t flags = _flags;
	uint32_t length;
	uint8_t* bases;
	CigarAlign* align = mapping._align;

	/*check the sequence strand*/
	if (mapping._strand) {
		flags |= SAM_FSR;
	}

	/*get sequence length*/
	length = _colorspace ? mapping._seqLength : seq._length;

	/*print out query name, bitwise-flag, and reference sequence name*/
	fprintf(_file, "%s\t%d\t%s\t", seq._name, flags,
			_genome->getGenomeName(mapping._genomeIndex));

	//print 1-based leftmost mapping position, mapping quality (phred-scale)*/
	fprintf(_file, "%u\t%d\t", (uint32_t) mapping._position, mapping._mapQual);

	//print extended CIGAR
	if (align->getCigar()) {
		align->cigarOut(_file);
	} else {
		fprintf(_file, "%dM", length);
	}
	fputc('\t', _file);

	//print paired-end information, INAVAILABLE
	fprintf(_file, "*\t0\t0\t");

	//print the query sequence
	if (_colorspace) {
		bases = (mapping._strand == 0) ? mapping._data : mapping._data + length;
	} else {
		bases = (mapping._strand == 0) ? seq._bases : seq._rbases;
	}

	for (uint32_t i = 0; i < length; ++i) {
		fputc(decode(bases[i]), _file);
	}
	fputc('\t', _file);

	//print the quality scores if available
	if (seq._quals) {
		uint8_t* quals;
		if (mapping._strand == 0) {
			quals = _colorspace ? mapping._data + 2 * length : seq._quals;
			for (uint32_t i = 0; i < length; ++i) {
				fputc(*quals, _file);
				++quals;
			}
		} else {
			quals = _colorspace ?
					mapping._data + 3 * length - 1 :
					seq._quals + seq._length - 1;
			for (uint32_t i = 0; i < length; ++i) {
				fputc(*quals, _file);
				--quals;
			}
		}
	} else {
		fputc('*', _file);
	}

	/*print tags*/
	fprintf(_file, "\tPG:Z:%s", PROGRAM_NAME);
	if (_groupID) {
		fprintf(_file, "\tRG:Z:%s", _groupID);
	}
	if(_groupLB){
		fprintf(_file, "\tLB:Z:%s", _groupLB);
	}
	fprintf(_file, "\tNM:i:%d", align->getEditDistance());
	fprintf(_file, "\tAS:i:%d", align->getAlignScore());

	/*end of line*/
	fputc('\n', _file);

}
/*print out the alignment with itself unaligned and its mate aligned*/
void SAM::printPaired(Sequence& seq, Mapping& mate, int _flags) {
	int flags = _flags;

	/*set the unmapped flag*/
	flags |= SAM_FSU;

	//check the strand of the mate
	if (mate._strand) {
		flags |= SAM_FMR;
	}

	/*print out query name, bitwise-flag, and reference sequence name*/
	fprintf(_file, "%s\t%d\t*\t", seq._name, flags);

	//print 1-based leftmost mapping position, mapping quality (phred-scale)*/
	fprintf(_file, "0\t0\t");

	//print extended CIGAR
	fputs("*\t", _file);

	/****************************
	 print the mate information
	 1-base mate mapping position and estimated insert size)
	 *****************************/
	//print the reference sequence mapped by the mate
	fprintf(_file, "%s\t", _genome->getGenomeName(mate._genomeIndex));

	//print the mapping position of the mate and the distance

	int64_t distance = 0; /*means not properly paired*/
	//Utils::log("self %d mate %d length %d strand %d-%d\n", self._position, mate._position, length, self._strand, mate._strand);

	fprintf(_file, "%d\t%ld\t", (int) mate._position, distance);

	//print the query sequence
	uint8_t* bases = seq._bases;
	for (uint32_t i = 0; i < seq._length; ++i) {
		fputc(decode(bases[i]), _file);
	}
	fputc('\t', _file);

	//print the quality scores if available
	if (seq._quals) {
		uint8_t* quals = seq._quals;
		for (uint32_t i = 0; i < seq._length; ++i) {
			fputc(*quals, _file);
			++quals;
		}
	} else {
		fputc('*', _file);
	}
	/*print tags*/
	fprintf(_file, "\tPG:Z:%s", PROGRAM_NAME);
	if (_groupID) {
		fprintf(_file, "\tRG:Z:%s", _groupID);
	}
	if(_groupLB){
		fprintf(_file, "\tLB:Z:%s", _groupLB);
	}
	fputc('\n', _file);
}
void SAM::printPaired(Sequence& seq, Mapping& self, Mapping& mate,
		bool properlyPaired, int _flags) {
	int flags = _flags;
	uint32_t length;
	uint8_t* bases;
	CigarAlign* align = self._align;

	/*the reads are paried*/
	if (properlyPaired) {
		flags |= SAM_FPP;
	}

	//check the strand of itself
	if (self._strand) {
		flags |= SAM_FSR;
	}

	//check the strand of the mate
	if (mate._strand) {
		flags |= SAM_FMR;
	}

	/*get sequence length*/
	length = _colorspace ? self._seqLength : seq._length;

	//print query-name, bitwise-flag, and reference-sequence-name
	fprintf(_file, "%s\t%d\t%s\t", seq._name, flags,
			_genome->getGenomeName(self._genomeIndex));

	//print 1-based leftmost mapping position, mapping quality (phred-scaled)
	fprintf(_file, "%d\t%d\t", (int) self._position, self._mapQual);

	//print extended CIGAR if applicable
	if (align->getCigar()) {
		align->cigarOut(_file);
	} else {
		fprintf(_file, "%dM", length);
	}
	fputc('\t', _file);

	/****************************
	 print the mate information
	 1-base mate mapping position and estimated insert size)
	 *****************************/
	//print the reference sequence mapped by the mate
	fprintf(_file, "%s\t",
			(self._genomeIndex == mate._genomeIndex) ?
					"=" : _genome->getGenomeName(mate._genomeIndex));

	//print the mapping position of the mate and the distance

	int64_t distance = 0; /*means not properly paired*/
	if (properlyPaired) {
		if (_pairingMode == MATE_PAIRED_READS) {
			distance = (int64_t) self._position - mate._position;
		} else {
			if (self._strand == 0) {
				distance = (int64_t) self._position - mate._position - length;
			} else {
				distance = (int64_t) self._position + length - mate._position;
			}
		}
	}
	//Utils::log("self %d mate %d length %d strand %d-%d\n", self._position, mate._position, length, self._strand, mate._strand);

	fprintf(_file, "%d\t%ld\t", (int) mate._position, distance);

	//print the query sequence on the same strand as the reference sequence
	if (_colorspace) {
		bases = (self._strand == 0) ? self._data : self._data + length;
	} else {
		bases = (self._strand == 0) ? seq._bases : seq._rbases;
	}
	for (uint32_t i = 0; i < length; ++i) {
		fputc(decode(bases[i]), _file);
	}
	fputc('\t', _file);

	//print the quality scores if available
	if (seq._quals) {
		uint8_t* quals;
		if (self._strand == 0) {
			quals = _colorspace ? self._data + 2 * length : seq._quals;
			for (uint32_t i = 0; i < length; ++i) {
				fputc(*quals, _file);
				++quals;
			}
		} else {
			quals = _colorspace ?
					self._data + 3 * length - 1 : seq._quals + seq._length - 1;
			for (uint32_t i = 0; i < length; ++i) {
				fputc(*quals, _file);
				--quals;
			}
		}
	} else {
		fputc('*', _file);
	}

	/*print tags*/
	fprintf(_file, "\tPG:Z:%s", PROGRAM_NAME);
	if (_groupID) {
		fprintf(_file, "\tRG:Z:%s", _groupID);
	}
	if(_groupLB){
		fprintf(_file, "\tLB:Z:%s", _groupLB);
	}
	fprintf(_file, "\tNM:i:%d", align->getEditDistance());
	fprintf(_file, "\tAS:i:%d", align->getAlignScore());
	fputc('\n', _file);
}

