#include "SuffixArray.h"
#include "Utils.h"

SuffixArray::SuffixArray(const char* saFileName) {
	FILE* file;
	uint32_t intVar;

	//open the suffix array file
	file = fopen(saFileName, "rb");
	if (!file) {
		Utils::exit("Faield to open suffix array file: %s", saFileName);
	}

	//get the $ symbol position
	fread(&intVar, sizeof(uint32_t), 1, file);

	//skip
	fseek(file, 4 * sizeof(uint32_t), SEEK_CUR);

	//get the suffix array interval
	fread(&_saFactor, sizeof(uint32_t), 1, file);

	//get the sequence length
	fread(&_seqLength, sizeof(uint32_t), 1, file);

	//get the suffix array data size
	_saSize = (_seqLength + _saFactor) / _saFactor;
	if (_saSize < 1)
		_saSize = 1;

	//read the suffix array data
	_saPtr = new uint32_t[_saSize];
	_saPtr[0] = (uint32_t) -1;
	fread(_saPtr + 1, sizeof(uint32_t), _saSize - 1, file);

	//close the file
	fclose(file);

	//report the suffix array memory size
	Utils::log("Suffix array memory size: %g MB\n",
			(_saSize * sizeof(uint32_t)) / 1024.0 / 1024.0);
}
SuffixArray::~SuffixArray() {
	if (_saPtr) {
		delete[] _saPtr;
	}
}

uint32_t SuffixArray::getFactor(const char* saFileName) {
	FILE* file;
	uint32_t intVar;

	//open the suffix array file
	file = fopen(saFileName, "rb");
	if (!file) {
		Utils::exit("Failed to open suffix array file: %s", saFileName);
	}

	//skip
	fseek(file, 5 * sizeof(uint32_t), SEEK_SET);

	//get the suffix array factor
	fread(&intVar, sizeof(uint32_t), 1, file);

	//close the file
	fclose(file);

	//return the suffix array factor
	return intVar;
}

