/*
 * SingleEnd.cpp
 *
 *  Created on: Jan 18, 2012
 *      Author: yongchao
 */

#include "SingleEnd.h"
#include "SeqFileParser.h"

SingleEnd::SingleEnd(Options* options, Genome * genome, SAM* sam) {
	_options = options;
	_genome = genome;
	_sam = sam;

	/*get the number of threads*/
	_numThreads = _options->getNumThreads();

	/*create parameters for threads*/
	_threads.resize(_numThreads);
	for (int tid = 0; tid < _numThreads; ++tid) {
		_threads[tid] = new ThreadParams(tid, _options, _sam,
				new MemEngine(_options, _genome, _sam), this);
	}
}

SingleEnd::~SingleEnd() {
	for (size_t i = 0; i < _threads.size(); ++i) {
		delete _threads[i];
	}
	_threads.clear();
}

void SingleEnd::execute() {
	SeqFileParser *parser;
	vector<pthread_t> threadIDs;

	/*initialize variables*/
	threadIDs.resize(_numThreads);

	/*for paired-end alignment*/
	vector<pair<string, int> > &inputs = _options->getInputFileList();
	for (size_t file = 0; file < inputs.size(); ++file) {
		parser = new SeqFileParser(_options, inputs[file].first.c_str(),
				_numThreads > 1, inputs[file].second);

		/*create threads*/
		for (size_t tid = 0; tid < _threads.size(); ++tid) {
			/*set file parser*/
			_threads[tid]->setFileParser(parser);

			/*create the thread*/
			if (pthread_create(&threadIDs[tid], NULL, _threadFunc,
					_threads[tid]) != 0) {
				Utils::exit("Thread creating failed\n");
			}
		}
		/*wait for the completion of all threads*/
		for (size_t tid = 0; tid < _threads.size(); ++tid) {
			pthread_join(threadIDs[tid], NULL);
		}
		/*release the file parser*/
		delete parser;
	}
	threadIDs.clear();

	/*report the alignment information*/
	long sum = 0, sum2 = 0;
	for (size_t tid = 0; tid < _threads.size(); ++tid) {
		sum += _threads[tid]->_numReads;
		sum2 += _threads[tid]->_numAligned;
	}
	Utils::log("#Reads aligned: %ld / %ld (%.2f%%)\n", sum2, sum,
			((double) (sum2)) / sum * 100);
}
void* SingleEnd::_threadFunc(void* arg) {

	Sequence seq;
	Mapping* mapping;
	int64_t numReads = 0, numAligned = 0;
	ThreadParams* params = (ThreadParams*) arg;
	Options *options = params->_options;
	SAM *sam = params->_sam;
	SeqFileParser* parser = params->_parser;
	MemEngine* engine = params->_engine;
	double stime = Utils::getSysTime();
	double etime;
	bool aligned;
	int multi = options->getMaxMultiAligns();
	int minMapQual = options->getMinMapQual();
	/*reserve space for the vector*/
	vector<Mapping*> mappings;
	mappings.reserve(multi);

	while (1) {
		/*read a sequence*/
		if (!parser->getSeq(seq)) {
			break;
		}
		/*invoke the engine*/
		engine->align(seq, mappings);

		/*output the alignment*/
		options->lock();
		if (mappings.size() > 0) {
			aligned = false;
			for (size_t i = 0; i < mappings.size(); ++i) {
				mapping = mappings[i];
				if (mapping->_mapQual >= minMapQual) {
					sam->print(seq, *mapping);
					aligned = true;
				} else {
					sam->print(seq); /*as unaligned*/
				}
			}
			if (aligned) {
				++numAligned;
			}
		} else { /*unaligned*/
			sam->print(seq);
		}
		options->unlock();

		/*release the mapping*/
		for (size_t i = 0; i < mappings.size(); ++i) {
			delete mappings[i];
		}
		mappings.clear();

		/*statistical information*/
		++numReads;
		if (numReads % 100000 == 0) {
			etime = Utils::getSysTime();
			Utils::log("processed %ld reads by thread %d in %.2f seconds\n",
					numReads, params->_tid, etime - stime);
		}
	}
	/*return the alignment results*/
	params->_numAligned += numAligned;
	params->_numReads += numReads;

	return NULL;
}
