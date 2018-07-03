/*
 * Thread.h
 *
 *  Created on: Jan 24, 2012
 *      Author: yongchao
 */

#ifndef THREAD_H_
#define THREAD_H_
#include "Macros.h"
#include "MemEngine.h"
#include "SeqFileParser.h"

struct ThreadParams
{
	ThreadParams(int tid, Options* options, SAM* sam, MemEngine *engine,
			void* aligner) {
		_tid = tid;
		_options = options;
		_sam = sam;
		_engine = engine;
		_parser = NULL;
		_aligner = aligner;

		/*for recording*/
		_numAligned = _numPaired = _numReads = 0;
	}
	~ThreadParams() {
		if (_engine) {
			delete _engine;
		}
	}
	inline void setFileParser(SeqFileParser* parser) {
		_parser = parser;
	}
	inline void setFileParser(SeqFileParser* parser1, SeqFileParser* parser2) {
		_parser1 = parser1;
		_parser2 = parser2;
	}

	int _tid;
	Options* _options;
	SAM* _sam;
	MemEngine* _engine;
	void* _aligner;
	SeqFileParser* _parser;
	SeqFileParser* _parser1;
	SeqFileParser* _parser2;
	int64_t _numAligned;
	int64_t _numPaired;
	int64_t _numReads;
};

#endif /* THREAD_H_ */
