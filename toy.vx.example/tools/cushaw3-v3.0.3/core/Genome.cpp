/*
 * Genome.cpp
 *
 *  Created on: Dec 28, 2011
 *      Author: yongchao
 */

#include "Genome.h"
#include "Utils.h"

void Genome::_init(string bwtFileName, string saFileName, string annFileName,
		string baseBitmapFileName, string pacFileName, string basePacFileName,
		bool colorspace, bool maskAmbiguous) {
	FILE* file;
	char buffer[4096];
	int intVar;
	long long longVar;
	BwtAnn* ann;
	size_t length;

	//read .bwt file
	_bwt = new BWT(bwtFileName.c_str());

	//read .sa file
	_sa = new SuffixArray(saFileName.c_str());

	//read .ann file
	file = fopen(annFileName.c_str(), "rb");
	if (file == NULL) {
		Utils::exit("Failed to open file %s\n", annFileName.c_str());
	}
	fscanf(file, "%lld%d%u", &longVar, &_numSeqs, &_seed);

	_genomeLength = longVar;
	_anns = new BwtAnn[_numSeqs];
	if (_anns == NULL) {
		Utils::exit("Memory allocation failed in function %s in line %d\n",
				__FUNCTION__, __LINE__);
	}
	for (int i = 0; i < _numSeqs; ++i) {
		//get the pointer address

		ann = _anns + i;
		//read gi and sequence name
		fscanf(file, "%u%s", &(ann->_gi), buffer);
		length = strlen(buffer);
		ann->_name = new char[length + 1];
		if (ann->_name == NULL) {
			Utils::exit("Memory allocation failed in function %s in line %d\n",
					__FUNCTION__, __LINE__);
		}
		strcpy(ann->_name, buffer);

		/*read comments*/
		int c;
		char* p = buffer;
		while ((c = fgetc(file)) != '\n' && c != EOF)
			*p++ = c;
		*p = '\0';

		if (p - buffer > 1) {
			length = strlen(buffer);
			ann->_anno = new char[length];
			if (ann->_anno == NULL) {
				Utils::exit(
						"Memory allocation failed in function %s in line %d\n",
						__FUNCTION__, __LINE__);
			}
			strcpy(ann->_anno, buffer + 1);
		} else {
			ann->_anno = new char[1];
			if (ann->_anno == NULL) {
				Utils::exit(
						"Memory allocation failed in function %s in line %d\n",
						__FUNCTION__, __LINE__);
			}
			ann->_anno[0] = '\0';
		}

		/*read the remaining part*/
		fscanf(file, "%lld%d%d", &longVar, &ann->_length, &ann->_numAmbs);
		ann->_offset = longVar;
	}
	fclose(file);

#if 0
	/*print out the sequence information in the genome*/
	for (int i = 0; i < _numSeqs; ++i)
	Utils::log("@SQ\tSN:%s\tLN:%d\n", _anns[i]._name,
			_anns[i]._length);
#endif

	/*read .pac file*/
	file = fopen(pacFileName.c_str(), "rb");
	if (file == NULL) {
		Utils::exit("Failed to open file %s\n", pacFileName.c_str());
	}

	/*calculate the paced genome size*/
	size_t pacGenomeBytes = (size_t) ((_genomeLength + 3) >> 2);

	/*allocate space*/
	_pacGenome = new uint8_t[pacGenomeBytes];
	if (_pacGenome == NULL) {
		Utils::exit("Memory allocation failed in function %s line %d\n",
				__FUNCTION__, __LINE__);
	}
	//read the file
	size_t bytes = fread(_pacGenome, 1, pacGenomeBytes, file);
	if (bytes != pacGenomeBytes) {
		Utils::exit("Incomplete file %s (read %ld != %ld )\n",
				pacFileName.c_str(), bytes, pacGenomeBytes);
	}
	fclose(file);

	/*color-space reads?*/
	if (colorspace) {
		/*read the base-space .pac file, if applicable*/
		file = fopen(basePacFileName.c_str(), "rb");
		if (file == NULL) {
			Utils::exit("Failed to open file %s\n", basePacFileName.c_str());
		}
		/*allocate space*/
		_basePacGenome = new uint8_t[pacGenomeBytes];
		if (_basePacGenome == NULL) {
			Utils::exit("Memory allocation failed in function %s line %d\n",
					__FUNCTION__, __LINE__);
		}
		//read the file
		bytes = fread(_basePacGenome, 1, pacGenomeBytes, file);
		if (bytes != pacGenomeBytes) {
			Utils::exit("Incomplete file %s (read %ld != %ld )\n",
					basePacFileName.c_str(), bytes, pacGenomeBytes);
		}
		fclose(file);
	}

	if (maskAmbiguous) {
		/*open the bitmap file*/
		file = fopen(baseBitmapFileName.c_str(), "rb");
		if (file == NULL) {
			Utils::exit("Failed to open file %s\n", baseBitmapFileName.c_str());
		}
		/*allocate space*/
		size_t baseBitmapBytes = (_genomeLength + 7) >> 3;
		_baseBitmap = new uint8_t[baseBitmapBytes];
		if (_baseBitmap == NULL) {
			Utils::exit("Memory allocation failed in function %s line %d\n",
					__FUNCTION__, __LINE__);
		}
		//read the file
		bytes = fread(_baseBitmap, 1, baseBitmapBytes, file);
		if (bytes != baseBitmapBytes) {
			Utils::exit("Incomplete file %s (read %ld != %ld )\n",
					baseBitmapFileName.c_str(), bytes, baseBitmapBytes);
		}
		fclose(file);
	}
}

void Genome::getGenomeIndex(uint32_t position, int& genomeIndex) {
	int left, mid, right;
	if (position >= _genomeLength) {
		Utils::exit(
				"getGenomeIndex: Mapping position is longer than sequence length (%lld >= %lld)",
				position, _genomeLength);
	}
	// binary search for the sequence ID. Note that this is a bit different from the following one...

	left = 0;
	mid = 0;
	right = _numSeqs;
	while (left < right) {
		mid = (left + right) >> 1;
		if (position >= _anns[mid]._offset) {
			if (mid == _numSeqs - 1)
				break;
			if (position < _anns[mid + 1]._offset)
				break;

			left = mid + 1;
		} else {
			right = mid;
		}
	}
	genomeIndex = mid;
}

