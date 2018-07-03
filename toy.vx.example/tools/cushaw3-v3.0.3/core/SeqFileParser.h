/*
 * SeqFileParser.h
 *
 *  Created on: Dec 23, 2011
 *      Author: yongchao
 */

#ifndef SEQFILEPARSER_H_
#define SEQFILEPARSER_H_

#include "Macros.h"
#include "Sequence.h"
#include "Utils.h"
#include "MyFile.h"
#include "Options.h"
extern "C" {
#include "sam.h"
}

class SeqFileParser
{
public:
	/*public member functions*/
	SeqFileParser(Options* options, const char* path, bool withLock, int type,
			size_t BUFFER_SIZE = 4095);
	~SeqFileParser();

	//get the next sequence from the file
	inline size_t getSeq(Sequence& seq) {
		size_t ret;
		int32_t numNs = 0, i;
		/*read the sequence from the file*/
		_lock();
		if (_format == FILE_FORMAT_FASTA) {
			ret = getFastaSeq(seq);
		} else if (_format == FILE_FORMAT_FASTQ) {
			ret = getFastqSeq(seq);
		} else { /*BAM/SAM format*/
			ret = getBSamSeq(seq);
		}
		_unlock();

		/*compute the reverse complement*/
		if (ret > 0) {
			/*trim all consecutive Ns at the 3 prime end for each read*/
			for (i = seq._length - 1; i >= 0; --i){
				if(seq._bases[i] != UNKNOWN_BASE){
					break;
				}
			}
			seq._length = max(1, i + 1);

			/*compute the number of Ns in the full length*/
			for (i = seq._length - 1; i >= 0; --i) {
				if (seq._bases[i] == UNKNOWN_BASE) {
					++numNs;
				}
			}
			seq._tlength = seq._length - numNs;

			//compute the reverse complement of the sequence
			reverseComp(seq._rbases, seq._bases, seq._length, _colorspace);
		}
		return ret;
	}
	inline size_t getSeqLockFree(Sequence& seq) {
		size_t ret;
		int32_t numNs = 0, i;

		/*read the sequence from the file*/
		if (_format == FILE_FORMAT_FASTA) {
			ret = getFastaSeq(seq);
		} else if (_format == FILE_FORMAT_FASTQ) {
			ret = getFastqSeq(seq);
		} else { /*BAM/SAM format*/
			ret = getBSamSeq(seq);
		}

		/*compute the reverse complement*/
		if (ret > 0) {
			/*trim all consecutive Ns at the 3 prime end for each read*/
			for (i = seq._length - 1; i >= 0; --i){
				if(seq._bases[i] != UNKNOWN_BASE){
					break;
				}
			}
			seq._length = max(1, i + 1);

			/*compute the number of Ns in the full length*/
			for (i = seq._length - 1; i >= 0; --i) {
				if (seq._bases[i] == UNKNOWN_BASE) {
					++numNs;
				}
			}
			seq._tlength = seq._length - numNs;

			//compute the reverse complement of the sequence
			reverseComp(seq._rbases, seq._bases, seq._length, _colorspace);
		}

		return ret;
	}

	inline size_t getSeqLockFree(Sequence* seqs, int maxSeqs) {
		int index;
		for (index = 0; index < maxSeqs; ++index) {
			if (!getSeqLockFree(seqs[index])) {
				break;
			}
		}

		return index;
	}

	static inline void encode(uint8_t* s, size_t length) {
		uint8_t ch;
		for (size_t i = 0; i < length; i++) {
			ch = s[i];
			if (ch >= 'A' && ch <= 'Z') {
				ch -= 'A';
			} else if (ch >= 'a' && ch <= 'z') {
				ch -= 'a';
			} else {
				Utils::exit("Unexpected character %c at line %d in file %s\n",
						ch, __LINE__, __FILE__);
			}
			s[i] = _codeTab[ch];
		}
	}

private:
	/*private member functions*/
	void resizeBuffer(size_t nsize);
	size_t getFastaSeq(Sequence& seq);
	size_t getFastqSeq(Sequence& seq);
	size_t getBSamSeq(Sequence & seq);
	inline void reverseComp(uint8_t* rbases, uint8_t* bases, size_t length,
			bool colorspace) {
		size_t off;
		size_t halfLength = length / 2;

		if (colorspace) {
			for (size_t i = 0; i < halfLength; i++) {
				off = length - i - 1;
				rbases[off] = bases[i];
				rbases[i] = bases[off];
			}
			if (length & 1) {
				rbases[halfLength] = bases[halfLength];
			}
		} else {
			for (size_t i = 0; i < halfLength; i++) {
				off = length - i - 1;
				rbases[off] = _complements[bases[i]];
				rbases[i] = _complements[bases[off]];
			}
			if (length & 1) {
				rbases[halfLength] = _complements[bases[halfLength]];
			}
		}
	}
	inline void reverseComp(uint8_t* bases, size_t length, bool colorspace) {
		uint8_t ch;
		size_t off;
		size_t halfLength = length / 2;

		if (colorspace) {
			for (size_t i = 0; i < halfLength; i++) {
				off = length - i - 1;
				ch = bases[off];
				bases[off] = bases[i];
				bases[i] = ch;
			}
		} else {
			for (size_t i = 0; i < halfLength; i++) {
				off = length - i - 1;
				ch = bases[off];
				bases[off] = _complements[bases[i]];
				bases[i] = _complements[ch];
			}
			if (length & 1) {
				bases[halfLength] = _complements[bases[halfLength]];
			}
		}
	}

	inline void _lock() {
		if (_withLock) {
			pthread_mutex_lock(&_mutex);
		}
	}
	inline void _unlock() {
		if (_withLock) {
			pthread_mutex_unlock(&_mutex);
		}
	}

	/*buffered file operations*/
	inline int myfgetc(MyFilePt file) {
		/*check the end-of-file*/
		if (_fileBufferSentinel >= _fileBufferLength) {
			/*re-fill the buffer*/
			_fileBufferSentinel = 0;
			/*read file*/
			_fileBufferLength = myfread(_fileBuffer, 1, 4096, file);
			if (_fileBufferLength == 0) {
				/*reach the end of the file*/
				if (myfeof(file)) {
					return -1;
				} else {
					Utils::exit("File reading failed in function %s line %d\n",
							__FUNCTION__, __LINE__);
				}
			}
		}
		/*return the current character, and increase the sentinel position*/
		return _fileBuffer[_fileBufferSentinel++];
	}
	inline int myungetc(int ch, MyFilePt file) {
		if (_fileBufferSentinel >= 0) {
			_fileBuffer[--_fileBufferSentinel] = ch;
		} else {
			Utils::log("Two consecutive ungetc operations occurred\n");
			return -1; /*an error occurred, return end-of-file marker*/
		}
		return ch;
	}
private:
	/*private member variables*/
	//buffer for file reading
	uint8_t* _buffer;
	size_t _length;
	size_t _size;

	//FASTA/FASTQ file handler
	MyFilePt _fp;
	uint8_t* _fileBufferR;
	uint8_t* _fileBuffer;
	int _fileBufferLength;
	int _fileBufferSentinel;

	//BAM/SAM file handler
	samfile_t* _samfp;
	bam1_t* _samentry;
	//
	int _format;
	bool _colorspace;
	bool _trimPrimer;

	/*interal lock for concurrent accesses to a single file*/
	pthread_mutex_t _mutex;
	bool _withLock;

	static const uint8_t _codeTab[26];
	static const uint8_t _decodeTab[5];
	static const uint8_t _complements[5];
};

#endif /* SEQFILEPARSER_H_ */
