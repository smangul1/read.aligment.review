/*
 * BSOptions.h
 *
 *  Created on: Mar 19, 2013
 *      Author: yongchao
 */

#ifndef BSOPTIONS_H_
#define BSOPTIONS_H_
#include "Options.h"
class BSOptions : public Options
{
public:
	BSOptions();
	~BSOptions(){}

	/*override functions*/
	void printUsage();
	bool parse(int argc, char* argv[]);

protected:
	void _setDefaults();
};




#endif /* BSOPTIONS_H_ */
