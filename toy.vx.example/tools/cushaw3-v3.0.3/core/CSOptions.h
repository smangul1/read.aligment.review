/*
 * CSOptions.h
 *
 *  Created on: Mar 19, 2013
 *      Author: yongchao
 */

#ifndef CSOPTIONS_H_
#define CSOPTIONS_H_
#include "Options.h"

class CSOptions : public Options
{
public:
	CSOptions();
	~CSOptions(){}

	/*override functions*/
	void printUsage();
	bool parse(int argc, char* argv[]);

protected:
	void _setDefaults();
};


#endif /* CSOPTIONS_H_ */
