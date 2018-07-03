/*
 * Genome.h
 *
 *  Created on: Dec 28, 2011
 *      Author: yongchao
 */

#ifndef GENOME_H_
#define GENOME_H_
#include "Macros.h"
#include "BWT.h"
#include "SuffixArray.h"
#include "Utils.h"
#include "Options.h"

struct BwtAnn
{
	BwtAnn() {
		_name = _anno = NULL;
	}
	~BwtAnn() {
		if (_name)
			delete[] _name;
		if (_anno)
			delete[] _anno;
	}

	int64_t _offset;
	int32_t _length;
	int32_t _numAmbs;
	uint32_t _gi;
	char* _name;
	char* _anno;
};

class Genome
{
public:
	Genome(Options* options) {
		string& bwtBase = options->getBwtFileBase();
		string bwtFileName, saFileName, annFileName;
		string pacFileName, basePacFileName;
		string baseBitmapFileName;

		/*initalize all pointers*/
		_bwt = NULL;
		_sa = NULL;
		_anns = NULL;
		_pacGenome = NULL;
		_basePacGenome = NULL;
		_baseBitmap = NULL;

		/*using the reverse orientation of the genome*/
		bwtFileName = bwtBase + ".rbwt";
		saFileName = bwtBase + ".rsa";
		annFileName = bwtBase + ".ann";
		baseBitmapFileName = bwtBase + ".map";
		pacFileName = bwtBase + ".pac";
		basePacFileName = bwtBase + ".nt.pac"; /*for color-space*/

		/*read the data from files*/
		_init(bwtFileName, saFileName, annFileName, baseBitmapFileName,
				pacFileName, basePacFileName, options->isColorSpace(),
				options->maskAmbiguous());
	}
	~Genome() {
		if (_bwt)
			delete _bwt;
		if (_sa)
			delete _sa;
		if (_anns)
			delete[] _anns;
		if (_pacGenome) {
			delete[] _pacGenome;
		}
		if (_basePacGenome) {
			delete[] _basePacGenome;
		}
		if (_baseBitmap) {
			delete[] _baseBitmap;
		}
	}

	inline BWT* getBWT() {
		return _bwt;
	}
	inline SuffixArray* getSuffixArray() {
		return _sa;
	}
	inline uint8_t* getBaseBitmap() {
		return _baseBitmap;
	}
	inline uint8_t* getPacGenome() {
		return _pacGenome;
	}
	/*only applicable for color-space reads*/
	inline uint8_t* getBasePacGenome() {
		return _basePacGenome;
	}
	inline void refineRegionRange(int genomeIndex, int64_t& lowerBound,
			int64_t& upperBound) {
		int64_t left = _anns[genomeIndex]._offset;
		int64_t right = left + _anns[genomeIndex]._length - 1;

		lowerBound = max(lowerBound, left);
		lowerBound = min(lowerBound, right);

		upperBound = max(upperBound, left);
		upperBound = min(upperBound, right);
	}
	/*get the genome sequence index in the packed genome*/
	void getGenomeIndex(uint32_t position, int& genomeIndex);

	/*check if the alignment lies in a random filed*/
	int getRandomField(uint32_t position, int length);

	/*get the offset of the genome*/
	inline int64_t getGenomeOffset(int genomeIndex) {
		return _anns[genomeIndex]._offset;
	}

	/*get the genome sequence name*/
	inline char* getGenomeName(int genomeIndex) {
		return _anns[genomeIndex]._name;
	}

	/*print out all genome names in SAM format*/
	inline void genomeNamesOut(FILE* file) {
		if (file == NULL) {
			Utils::exit("Invalid file pointer in function %s line %d\n",
					__FUNCTION__, __LINE__);
		}
		for (int i = 0; i < _numSeqs; ++i) {
			fprintf(file, "@SQ\tSN:%s\tLN:%d\n", _anns[i]._name,
					_anns[i]._length);
		}
	}
	inline void genomeNamesOut(gzFile file) {
		if (file == NULL) {
			Utils::exit("Invalid file pointer in function %s line %d\n",
					__FUNCTION__, __LINE__);
		}
		for (int i = 0; i < _numSeqs; ++i) {
			gzprintf(file, "@SQ\tSN:%s\tLN:%d\n", _anns[i]._name,
					_anns[i]._length);
		}
	}

private:
	/*private member variables*/
	BWT* _bwt;
	SuffixArray* _sa;

	/*genome information*/
	int64_t _genomeLength; /*the original genome length*/
	uint8_t* _pacGenome;
	uint8_t* _basePacGenome; /*if it is color space, this genome is for base-space; otherwise, NULL*/
	uint8_t* _baseBitmap; /*bitmap for unknown bases*/

	uint32_t _seed; /*the fixed random number generator seed*/
	BwtAnn* _anns; // _numSeqs elements
	int _numSeqs;

private:
	/*private member functions*/
	void _init(string bwtFileName, string saFileName, string annFileName,
			string baseBitmapFileName, string pacFileName,
			string basePacFileName, bool colorspace, bool maskAmbiguous);
};
#endif /* GENOME_H_ */
