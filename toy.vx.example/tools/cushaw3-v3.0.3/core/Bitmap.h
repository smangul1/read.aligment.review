/*
 * Bitmap.h
 *
 *  Created on: Jan 6, 2013
 *      Author: yongchao
 */

#ifndef BITMAP_H_
#define BITMAP_H_
#include "Macros.h"
#include "Utils.h"

class Bitmap
{
public:
	Bitmap(size_t size) {
		_dataSize = max(size, (size_t)1);
		_data = new uint8_t[size];
		if (!_data) {
			Utils::exit("Memory allocation failed at line %d in function %s\n",
					__LINE__, __FUNCTION__);
		}

		/*clear all bits*/
		memset(_data, 0, _dataSize);

		_numBitsSet = 0;
	}
	~Bitmap() {
		if (_data) {
			delete[] _data;
		}
	}

	inline int64_t getNumBitsSet()
	{
		return _numBitsSet;
	}
	inline void set() {
		if (_data) {
			memset(_data, 0x0FF, _dataSize);
		}
	}
	inline void set(size_t position) {

		/*resize the buffer*/
		if (position >= _dataSize) {
			resize(_dataSize);
		}
		/*set the position*/
		size_t bytes = position >> 3;
		size_t bits = position & 7;
		uint8_t data = _data[bytes];

		/*if this bit has not been set*/
		if(!((data >> bits) & 1)){
			++_numBitsSet;
			_data[bytes] = data | (1 << bits);
		}
	}
	inline void reset() {
		if (_data) {
			memset(_data, 0, _dataSize);
		}
		_numBitsSet = 0;
	}
	inline void reset(size_t position) {
		/*resize the buffer*/
		if (position >= _dataSize) {
			resize(_dataSize / 2);
		}

		/*set the position*/
		size_t bytes = position >> 3;
		size_t bits = position & 7;
		uint8_t data = _data[bytes];
		if((data >> bits) & 1){
			--_numBitsSet;
			_data[bytes] = data & (~(1 << bits));
		}
	}
	inline bool test(size_t position) {
		/*resize the buffer*/
		if (position >= _dataSize) {
			Utils::exit("Positions (%ld) exceeds the bitmap size (%ld)\n", position, _dataSize);
		}
		/*set the position*/
		size_t bytes = position >> 3;
		size_t bits = position & 7;

		return (_data[bytes] & (1 << bits));
	}
private:
	inline void resize(size_t inc) {
		if (inc == 0) {
			return;
		}

		if (_data) {
			/*re-allocate memory*/
			uint8_t* buffer = new uint8_t[_dataSize + inc];
			if (!buffer) {
				Utils::exit(
						"Memory allocation failed at line %d in function %s\n",
						__LINE__, __FUNCTION__);
			}
			memcpy(buffer, _data, _dataSize);
			memset(buffer + _dataSize, 0, inc);

			delete[] _data;
			_data = buffer;
			_dataSize += inc;
		} else {
			/*new buffer*/
			_dataSize = inc;
			_data = new uint8_t[_dataSize];
			if (!_data) {
				Utils::exit(
						"Memory allocation failed at line %d in function %s\n",
						__LINE__, __FUNCTION__);
			}
			memset(_data, 0, _dataSize);
		}
	}
	uint8_t* _data;
	size_t _dataSize;
	int64_t _numBitsSet;
};

#endif /* BITMAP_H_ */
