/*
 * Options.h
 *
 *  Created on: Jan 10, 2012
 *      Author: yongchao
 */

#ifndef OPTIONS_H_
#define OPTIONS_H_
#include "Macros.h"
#include "Utils.h"

class Options
{
public:
	Options();
	virtual ~Options();

	/*virtual functions*/
	virtual void printUsage() = 0;
	virtual bool parse(int argc, char* argv[]) = 0;

	//get the file name of the suffix array index
	inline string& getSamFileName() {
		return _samFileName;
	}
	//get the prefix of BWT file name
	inline string& getBwtFileBase() {
		return _bwtFileBase;
	}
	//get the input file list
	inline vector<pair<string, int> >& getInputFileList() {
		return _readsFileNames;
	}
	inline int getMaxMultiAligns() {
		return _maxMultiAligns;
	}
	inline bool isSensitiveParing()
	{
		return _sensitivePairing;
	}
	inline bool isPaired() {
		return _paired;
	}
	inline bool rescueMate() {
		return _rescueMate;
	}
	inline int getPairingMode()
	{
		return _pairingMode;
	}
	inline bool isColorSpace() {
		return _colorspace;
	}
	inline bool trimPrimer()
	{
		return _trimPrimer;
	}
	inline bool maskAmbiguous()
	{
		return _maskAmbiguous;
	}
	inline int rescueTwice()
	{
		return _rescueTwice;
	}
	inline void setPaired(bool pe) {
		_paired = pe;
	}
	inline uint32_t getMaxSeedOcc() {
		return _maxSeedOcc;
	}
	inline float getMinRatio() {
		return _minRatio;
	}
	inline float getMinIdentity() {
		return _minIdentity;
	}
	inline float getGMinIdentity() {
		return _gminIdentity;
	}
	inline int getMinAlignScore() {
		return _minAlginScore;
	}
	inline bool estimateInsertSize() {
		return _estInsertSize;
	}
	inline int getInsertSize() {
		return _insertSize;
	}
	inline void setInsertSize(int insertSize) {
		_insertSize = insertSize;
	}
	inline int getStdInsertSize() {
		return _stdInsertsize;
	}
	inline void setStdInsertSize(int stdInsertSize) {
		_stdInsertsize = stdInsertSize;
	}
	inline void setTopReadsEstIns(int nreads) {
		_topReadsEstIns = nreads;
	}
	inline int getTopReadsEstIns() {
		return _topReadsEstIns;
	}
	inline void setMapQualEstIns(int mapQ) {
		_mapQualReliable = mapQ;
	}
	inline int getMapQualReliable() {
		return _mapQualReliable;
	}
	inline int getMinMapQual()
	{
		return _minMapQual;
	}
	inline int getNumThreads() {
		return _numThreads;
	}
	inline int getGapOpen() {
		return _gapOpen;
	}
	inline int getGapExtend() {
		return _gapExtend;
	}
	inline int getMismatch() {
		return _mismatch;
	}
	inline int getMatch() {
		return _match;
	}
	inline int getAlignType() {
		return _alignType;
	}
	inline int getMaxSeedsPerBatch() {
		return _maxSeedsPerBatch;
	}
	inline int getMaxEditDistance()
	{
		return _maxEditDistance;
	}
	inline int getNumErrors(uint32_t length) {
		if (length > MAX_SEQ_LENGTH) {
			return estimateNumErrors(length, _missProb);
		}
		return _numErrorTable[length];
	}
	inline int getMaxSeedPairs() {
		return _maxSeedPairs;
	}
	inline int getMinSeedSize(uint32_t length) {
		if (length > MAX_SEQ_LENGTH) {
			int numErrors = estimateNumErrors(length, _missProb);
			_numErrorTable[length] = numErrors;

			/*estimate the minimal seed size according to dove hole principle*/
			int seedSize = length / (numErrors + 1);
			if (seedSize < _lowerMinSeedSize) {
				seedSize = _lowerMinSeedSize;
			}
			if (seedSize > _upperMinSeedSize) {
				seedSize = _upperMinSeedSize;
			}
			return seedSize;
		}
		return _minSeedSizeTable[length];
	}
	inline void lock() {
		if (_numThreads > 1) {
			pthread_mutex_lock(&globalMutex);
		}
	}
	inline void unlock() {
		if (_numThreads > 1) {
			pthread_mutex_unlock(&globalMutex);
		}
	}

protected:
	/*member variables*/
	string _bwtFileBase;
	string _samFileName;
	vector<pair<string, int> > _readsFileNames;
	uint32_t _maxSeedOcc;
	float _minRatio;
	float _minIdentity;
	float _gminIdentity;
	float _missProb;
	int _insertSize;
	int _stdInsertsize;
	int _topReadsEstIns;
	int _mapQualReliable;
	int _minMapQual;
	int _maxMultiAligns;
	int _numThreads;
	int _lowerMinSeedSize;
	int _upperMinSeedSize;
	int _maxSeedPairs;
	int _maxSeedsPerBatch;
	int _maxEditDistance;

	/*scoring scheme*/
	int _gapOpen;
	int _gapExtend;
	int _mismatch;
	int _match;
	int _minAlginScore;
	int _alignType;
	int _pairingMode;
	int _rescueTwice;

	int *_numErrorTable;
	int *_minSeedSizeTable;

	bool _estInsertSize;
	bool _sensitivePairing;
	bool _paired;
	bool _maskAmbiguous;
	bool _rescueMate;
	bool _colorspace;
	bool _trimPrimer;	/*trim the primer base*/

	/*read group information*/
	string _rgID; /*read group identifier*/
	string _rgSM; /*sample name. Use pool name where a pool is being sequence*/
	string _rgLB; /*library*/
	string _rgPL; /*platform/technology used to produce the reads*/
	string _rgPU; /*platform unit*/
	string _rgCN; /*name of sequencing center produced the read*/
	string _rgDS; /*description*/
	string _rgDT; /*date the run was produced*/
	string _rgPI; /*predicted median insert size*/

	/*global file lock*/
	pthread_mutex_t globalMutex;

protected:
	/*virtual function*/
	virtual void _setDefaults();

	/*estimate number of errors from read length and mean error rate*/
	inline int estimateNumErrors(size_t length, float errorRate) {
		float uniformErrorRate = 0.02;
		double elambda = exp(-length * uniformErrorRate);
		double sum, y = 1.0;
		int k, x = 1;
		for (k = 1, sum = elambda; k < 1000; ++k) {
			y *= length * uniformErrorRate;
			x *= k;
			sum += elambda * y / x;
			if (1.0 - sum < errorRate) {
				return k;
			}
		}
		return 2;
	}

	static int _atypes[2];
	/*set friend class*/
	friend class SAM;
};

#endif /* OPTIONS_H_ */
