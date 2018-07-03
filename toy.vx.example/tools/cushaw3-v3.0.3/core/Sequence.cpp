/*
 * Sequence.cpp
 *
 *  Created on: Dec 23, 2011
 *      Author: yongchao
 */

#include "Sequence.h"
#include "Utils.h"

Sequence::Sequence() {
	_name = NULL;
	_bases = NULL;
	_rbases = NULL;
	_quals = NULL;
	_nameSize = _seqSize = 0;
	_length = 0;
	_tlength = 0;
}
Sequence::Sequence(const Sequence & s) {
	_length = s._length;
	_tlength = s._tlength;
	_nameSize = s._nameSize;
	_seqSize = s._seqSize;
	if (_length == 0) {
		_name = NULL;
		_bases = NULL;
		_rbases = NULL;
		_quals = NULL;
		_nameSize = 0;
		_seqSize = 0;
		return;
	}
	if (s._name) {
		_name = new uint8_t[_nameSize];
		if (_name == NULL) {
			Utils::exit("Memory allocation failed in function %s line %d\n",
					__FUNCTION__, __LINE__);
		}
		strcpy((char*) _name, (const char*) s._name);
	}
	if (s._bases) {
		_bases = new uint8_t[_seqSize];
		if (_bases == NULL) {
			Utils::exit("Memory allocation failed in function %s line %d\n",
					__FUNCTION__, __LINE__);
		}
		memcpy(_bases, s._bases, _length);
	}
	if (s._rbases) {
		_rbases = new uint8_t[_seqSize];
		if (_rbases == NULL) {
			Utils::exit("Memory allocation failed in function %s line %d\n",
					__FUNCTION__, __LINE__);
		}
		memcpy(_rbases, s._rbases, _length);
	}
	if (s._quals) {
		_quals = new uint8_t[_seqSize];
		if (_quals == NULL) {
			Utils::exit("Memory allocation failed in function %s line %d\n",
					__FUNCTION__, __LINE__);
		}
		memcpy(_quals, s._quals, _length);
	}
}
Sequence::~Sequence() {
	clear();
}
void Sequence::clear() {
	if (_name) {
		delete[] _name;
	}
	if (_bases) {
		delete[] _bases;
	}
	if (_rbases) {
		delete[] _rbases;
	}
	if (_quals) {
		delete[] _quals;
	}

	_name = NULL;
	_bases = NULL;
	_rbases = NULL;
	_quals = NULL;
	_length = 0;
	_tlength = 0;
	_nameSize = 0;
	_seqSize = 0;
}
void Sequence::setNameSize(size_t size) {
	if (size >= _nameSize) {
		_nameSize = size * 2;
		if (_name) {
			delete[] _name;
		}
		_name = new uint8_t[_nameSize];
		if (_name == NULL) {
			Utils::exit("Memory allocation failed in function %s line %d\n",
					__FUNCTION__, __LINE__);
		}
	}
}
void Sequence::setSequenceSize(size_t size, bool quals) {
	if (size >= _seqSize) {
		_seqSize = size * 2;
		/*forward strand*/
		if (_bases) {
			delete[] _bases;
		}
		_bases = new uint8_t[_seqSize];
		if (_bases == NULL) {
			Utils::exit("Memory allocation failed in function %s line %d\n",
					__FUNCTION__, __LINE__);
		}
		/*reverse strand*/
		if (_rbases) {
			delete[] _rbases;
		}
		_rbases = new uint8_t[_seqSize];
		if (_rbases == NULL) {
			Utils::exit("Memory allocation failed in function %s line %d\n",
					__FUNCTION__, __LINE__);
		}

		/*allocate space for quality scores*/
		if (quals) {
			if (_quals) {
				delete[] _quals;
			}
			_quals = new uint8_t[_seqSize];
			if (_quals == NULL) {
				Utils::exit("Memory allocation failed in function %s line %d\n",
						__FUNCTION__, __LINE__);
			}
		}
	}
}
void Sequence::print(FILE* file) {
	//print the sequence name
	if (_quals) {
		fputc('@', file);
	} else {
		fputc('>', file);
	}
	fprintf(file, "%s\n", _name);

	//print the query sequence
	for (uint32_t i = 0; i < _length; ++i) {
		fputc(decode(_bases[i]), file);
	}
	fputc('\n', file);

	//print the quality scores if available
	if (_quals) {
		/*printout comments*/
		fputc('+', file);
		fputc('\n', file);
		for (uint32_t i = 0; i < _length; ++i) {
			fputc(_quals[i], file);
		}
	}
	fputc('\n', file);
}
