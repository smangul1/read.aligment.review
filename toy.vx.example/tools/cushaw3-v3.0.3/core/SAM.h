/*
 * SAM.h
 *
 *  Created on: Jan 5, 2012
 *      Author: yongchao
 */

#ifndef SAM_H_
#define SAM_H_
#include "Macros.h"
#include "Genome.h"
#include "CigarAlign.h"
#include "Sequence.h"
#include <pthread.h>
#include "Options.h"
#include "Mapping.h"

class SAM
{
public:
	SAM(Options* options, Genome* genome);
	~SAM();

	/*open the output file*/
	void open();
	/*close the output file*/
	void close();
	void print(Sequence& seq, int flags = 0);
	void print(Sequence& seq, Mapping& mapping, int flags = 0);
	void printPaired(Sequence& seq, Mapping& mate, int flags = 0);
	void printPaired(Sequence& seq, Mapping& self, Mapping& mate, bool properlyPaired, int _flags);

private:
	/*private member variables*/
	int _numThreads;
	string _fileName;
	Genome* _genome;
	bool _colorspace;
	bool _pairingMode;
	char* _groupID;
	char* _groupLB;

	FILE* _file;
};

#endif /* SAM_H_ */
