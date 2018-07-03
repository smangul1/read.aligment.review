/*
 * SingleEnd.h
 *
 *  Created on: Jan 18, 2012
 *      Author: yongchao
 */

#ifndef SINGLEEND_H_
#define SINGLEEND_H_
#include "Macros.h"
#include "Options.h"
#include "MemEngine.h"
#include "Thread.h"

class SingleEnd
{
public:
	SingleEnd(Options* options, Genome * genome, SAM* sam);
	~SingleEnd();

	/*execute the paired-end alignment*/
	void execute();

	/*get thread parameter*/
	inline ThreadParams* getThreadParams(int tid) {
		return _threads[tid];
	}
private:
	/*private member variables*/
	Options* _options;
	Genome* _genome;
	SAM* _sam;
	int _numThreads;
	/*thread parameters*/
	vector<ThreadParams*> _threads;
private:
	/*private member functions*/
	static void* _threadFunc(void*);
};

#endif /* SINGLEEND_H_ */
