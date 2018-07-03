/*
 * Options.cpp
 *
 *  Created on: Jan 10, 2012
 *      Author: yongchao
 */

#include "Options.h"
#include "Utils.h"

int Options::_atypes[]=
{
		HAVE_LOCAL_ALIGN,
		HAVE_LOCAL_SEMI_GLOBAL_ALIGN
};
Options::Options() {
}
Options::~Options() {
	_readsFileNames.clear();
	delete[] _numErrorTable;
	delete[] _minSeedSizeTable;
	pthread_mutex_destroy(&globalMutex);
}
void Options::_setDefaults() {
}

void Options::printUsage()
{

}
bool Options::parse(int argc, char* argv[])
{
	return true;
}
