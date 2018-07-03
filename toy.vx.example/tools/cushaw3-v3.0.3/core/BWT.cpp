#include "BWT.h"
#include "Utils.h"

BWT::BWT(const char* bwtFile) {
	_bwtSeqLength = 0;
	_bwtDollar = 0;
	_bwtSize = 0;
	_bwtPtr = NULL;

	//initialize the counting array
	for (int i = 0; i < BWT_NUM_OCC; i++) {
		_bwtCCounts[i] = 0;
	}

	//load the BWT data from files
	this->loadFromFile(bwtFile);
}
BWT::~BWT() {
	if (_bwtPtr) {
		delete[] _bwtPtr;
	}
}
void BWT::loadFromFile(const char* bwtFileName) {
	FILE* bwtFile;

	bwtFile = fopen(bwtFileName, "rb");
	if (bwtFile == NULL) {
		Utils::exit("Failed to open BWT file: %s", bwtFileName);
	}

	//get the bwt size
	fseek(bwtFile, 0, SEEK_END);
	_bwtSize = (ftell(bwtFile) - sizeof(uint32_t) * 5) >> 2;

	//allocate spaces for the BWT
	_bwtPtr = new uint32_t[_bwtSize];

	//get the position of $ symbol in BWT
	fseek(bwtFile, 0, SEEK_SET);
	fread(&_bwtDollar, sizeof(uint32_t), 1, bwtFile);

	//get the first occurrence of a nucleotide 
	fread(_bwtCCounts + 1, sizeof(uint32_t), BWT_NUM_NUCLEOTIDE, bwtFile);

	//get the sequence length
	_bwtSeqLength = _bwtCCounts[BWT_NUM_NUCLEOTIDE];

	//read the BWT
	fread(_bwtPtr, sizeof(uint32_t), _bwtSize, bwtFile);

	//close the file
	fclose(bwtFile);

#if 0
	Utils::log("BWT file name: %s\n", bwtFileName);
	Utils::log("BWT size: %u Bytes\n", _bwtSize * sizeof(uint32_t));
	Utils::log("BWT dollar position: %u\n", _bwtDollar);
	Utils::log("Sequence length: %u\n", _bwtSeqLength);
	Utils::log("CCounts: ");
	for(int i = 0; i <= BWT_NUM_NUCLEOTIDE; i++)
	{
		Utils::log("%u, ", _bwtCCounts[i]);
	}
	Utils::log("\n");
#endif
}
uint32_t BWT::bwtOcc(uint8_t base, uint32_t pos) {
	uint32_t n;
	uint32_t* ptr;

	if (pos == (uint32_t) -1)
		return 0;
	if (pos == _bwtSeqLength)
		return _bwtCCounts[base + 1] - _bwtCCounts[base];

	//for bases indexed greater equal than the $ symbol, the index is decreased by 1 because $ is removed from the final BWT
	if (pos >= _bwtDollar)
		--pos;

	//get the address of the interleaved cumulative occurrence
	ptr = bwtInterleavedOcc(pos, base);
	n = ptr[base]; //get the cumulative occurrence of the base

	//skip the addresses of the interleaved cumulative occurrences for the four bases
	ptr += BWT_NUM_NUCLEOTIDE;

	// calculate Occ up to the last pos/32
	uint32_t _pos = (pos >> 5) << 5;
	for (uint32_t ipos = (pos >> BWT_OCC_INTERVAL_SHIFT) * BWT_OCC_INTERVAL;
			ipos < _pos; ipos += 32, ptr += 2) {
		n += bwtOccAux((uint64_t) ptr[0] << 32 | ptr[1], base);
	}

	// calculate Occ
	n += bwtOccAux(
			((uint64_t) ptr[0] << 32 | ptr[1])
					& ~((1ull << ((~pos & 31) << 1)) - 1), base);
	if (base == 0)
		n -= ~pos & 31; // corrected for the masked bits

	return n;
}
;
