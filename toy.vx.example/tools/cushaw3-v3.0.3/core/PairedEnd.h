/*
 * PairedEnd.h
 *
 *  Created on: Jan 18, 2012
 *      Author: yongchao
 */

#ifndef PAIREDEND_H_
#define PAIREDEND_H_
#include "Macros.h"
#include "Options.h"
#include "MemEngine.h"
#include "Thread.h"

class PairedEnd
{
public:
	PairedEnd(Options* options, Genome * genome, SAM* sam);
	~PairedEnd();

	/*execute the paired-end alignment*/
	void execute();
	void estimateInsertSize(int minAlignedPairs, int maxReadBatchSize,
			int numReadBatchs);
	inline void lock() {
		if (_numThreads > 1) {
			pthread_mutex_lock(&_mutex);
		}
	}
	inline void unlock() {
		if (_numThreads > 1) {
			pthread_mutex_unlock(&_mutex);
		}
	}
private:
	/*private member variables*/
	Options* _options;
	Genome* _genome;
	SAM* _sam;
	int _numThreads;
	/*thread parameters*/
	vector<ThreadParams*> _threads;

	pthread_mutex_t _mutex;
private:
	static void* _threadFunc(void*);
};

#endif /* PAIREDEND_H_ */
