/*
 * MemEngine.h
 *
 *  Created on: Dec 24, 2011
 *      Author: yongchao
 */

#ifndef MEMENGINE_H_
#define MEMENGINE_H_
#include "Macros.h"
#include "Utils.h"
#include "Genome.h"
#include "Sequence.h"
#include "Aligner.h"
#include "SAM.h"
#include "Options.h"
#include "Mapping.h"
#include "Structs.h"

/*engine for maximal exact match (MEM) seed generation*/
class MemEngine
{
public:
	MemEngine(Options* options, Genome* rgenome, SAM* sam);
	~MemEngine();

	/*for single-end alignment*/
	bool align(Sequence& seq, Mapping*& bestMapping, int seedType =
			MAXIMAL_EXACT_MATCH_SEED);

	/*for multiple single-end alignment*/
	bool align(Sequence& seq, vector<Mapping*>& mappings, int seedType =
			MAXIMAL_EXACT_MATCH_SEED);

	/*multiple paired-end alignment */
	bool align(Sequence& seq1, Sequence& seq2, vector<Mapping*>& mappings1,
			vector<Mapping*>& mappings2);

	/*update distance information*/
	inline void updateDistance() {
		_maxDistance = _options->getInsertSize()
				+ 4 * _options->getStdInsertSize();
	}

protected:
	/*protected member functions*/
	/*compute the alignment for a sequence*/
	inline Mapping* _getAlignment(Aligner* aligner, Sequence& seq,
			size_t numSeeds, Seed* seeds) {
		int32_t bestHit;
		int mapQual;

		/*get the single-end alignment*/
		if (numSeeds > 0) {
			bestHit = _getBestHit(_aligner, seq, numSeeds, seeds, mapQual);

			/*get the alignment*/
			if (bestHit >= 0) {
				return _getAlignment(_aligner, seq, seeds[bestHit], mapQual);
			}
		}
		return NULL;
	}
	/*calculate the mate read mapping region*/
	inline void _getMateRegion(int64_t& genomeStart, int64_t& genomeEnd,
			int64_t selfMapPosition, int selfStrand, int selfLength,
			int mateLength) {
#if 0
		if (selfStrand == 0) {
			genomeStart = selfMapPosition + 1;
			genomeEnd = selfMapPosition + _maxDistance + mateLength;
		} else {
			genomeStart = selfMapPosition + selfLength - _maxDistance;
			genomeEnd = selfMapPosition + selfLength - 1;
		}
#else
		/*We found that the coordinate of the left read is not always smaller than the coordinate of the right end*/
		genomeStart = selfMapPosition + selfLength - _maxDistance;
		genomeEnd = selfMapPosition + _maxDistance + mateLength;
#endif
	}

	/*normal single-end alignment*/
	Mapping* _getAlignment(Aligner* aligner, Sequence& seq, Seed& seeds,
			int mapQual);

	Mapping* _getAlignment(Aligner* aligner, Sequence& seq, int strand,
			int window, int genomeIndex, uint32_t genomeStart,
			size_t genomeLength, int mapQual, bool rescue);
	Mapping* _getAlignment(Aligner* aligner, int alignType, uint8_t* genome,
			uint8_t* seqBases, uint32_t seqLength, int strand,
			int64_t genomeIndex, int64_t targetPosition, CigarAlign* align,
			int mapQual, float minIdentity);

	/*get the best seed hits*/
	int32_t _getBestHits(Aligner* aligner, Sequence& seq, size_t numSeeds,
			Seed* seeds, vector<int32_t>& bestHits, size_t& numTopSeeds,
			int& mapQual);

	/*get only one best hit*/
	int32_t _getBestHit(Aligner* aligner, Sequence& seq, size_t numSeeds,
			Seed* seeds, int& mapQual);

	/*compute the alignment scores*/
	size_t _getAlignmentScores(Aligner* aligner, Sequence& seq, size_t numSeeds,
			Seed* seeds);

	/*generate seeds for a single strand. Returns true if it has an exact match in the full length, and
	 * false otherwise*/
	uint32_t _genMEMSeeds(Sequence& seq, int strand, uint4* ranges,
			uint32_t minSeedSize);

	/*generate fixed-length striped seeds*/
	uint32_t _genStripedSeeds(Sequence& seq, int strand, uint4* ranges,
			uint32_t minSeedSize, uint32_t stride);

	/*calcualte suffix array intervals*/
	uint32_t _locateSAI(Sequence& seq, vector<uint4>& ranges,
			vector<uint4>& rranges, uint32_t& numRanges, uint32_t& numRranges,
			int seedType);

	/*locate the mapping positions of each seeds and perform voting on conditions*/
	void _locateSeeds(Sequence& seq, vector<Seed>& seeds, int seedType);

	/*get the quatified seeds*/
	inline void _getQualifiedSeeds(vector<Seed>&qualSeeds, vector<Seed>& seeds,
			vector<int32_t>& bestHits, bool append) {
		size_t oldSize;
		if (append == false) {
			qualSeeds.clear();
			qualSeeds.resize(bestHits.size());
			for (size_t i = 0; i < bestHits.size(); ++i) {
				qualSeeds[i] = seeds[bestHits[i]];
			}
		} else {
			/*resize memory*/
			oldSize = qualSeeds.size();
			qualSeeds.resize(oldSize + bestHits.size());
			for (size_t i = 0; i < bestHits.size(); ++i) {
				qualSeeds[oldSize + i] = seeds[bestHits[i]];
			}

			/*re-sort all seeds*/
			stable_sort(qualSeeds.begin(), qualSeeds.end());
		}
#if 0
		/*print out the alignment scores*/
		Utils::log("-------------Top qualified seeds:----------------\n");
		for (size_t i = 0; i < qualSeeds.size(); ++i) {
			Seed* p = &qualSeeds[i];
			fprintf(stderr, "%u %u %u %u %u (genome)%d\n", p->_queryPosition,
					p->_targetPosition, p->_seedLength, p->_strand,
					p->_alignScore, p->_genomeIndex);
		}
#endif
	}
	/*color space to base space conversion using BWA approach*/
	pair<uint8_t*, int32_t> _cs2bsReadConvertBWA(const Sequence& seq,
			Mapping* mapping);

	/*compute priority scores*/
	inline float _getHeapPriority(uint32_t leftLength, uint32_t rightLength,
			float leftScore, float rightScore) {
		float score, score2;

		score = leftScore < 1 ? 1 : leftScore;
		score2 = rightScore < 1 ? 1 : rightScore;
		score = score / (float) (leftLength * _matchScore);
		score2 = score2 / (float) (rightLength * _matchScore);
		//Utils::log("align scores: %g %g\n", score, score2);
		score = 2 * score * score2 / (score + score2);
		//score *= max(leftScore, rightScore);

		return score;
	}

	/*Color-space and base-space conversion*/
	inline uint8_t _genCSBase(uint8_t* bases, uint8_t* quals, int32_t strand,
			int32_t length, int32_t index) {
		int qual = quals[strand ? length - index - 1 : index] - 33;
		if (qual > 60) {
			qual = 60;
		}
		if (bases[index] > 3)
			qual = 63;
		return bases[index] << 6 | qual;
	}
	inline uint8_t _getRefBase(uint8_t* genome, int64_t gPosition) {
		return (genome[gPosition >> 2] >> ((~gPosition & 3) << 1)) & 3;
	}

	void _convertCS2BS(uint8_t* csRead, uint8_t* bsRead, uint8_t* bsRef,
			uint8_t* btArray, int32_t size);
	uint8_t* _genBSQuals(uint8_t *bsRead, uint8_t *csRead, uint8_t *tmpArray,
			int size);
private:
	/*private member variables*/
	Options* _options;
	Genome* _rgenome;
	BWT* _rbwt; /*pointer to the reverse BWT data of the target sequence*/
	SuffixArray* _rsa; /*pointer to the reverse Suffix Array data of the target sequence*/
	int64_t _bwtSeqLength;
	uint8_t* _baseBitmap;
	uint8_t* _pacGenome;
	uint8_t* _basePacGenome;

	/*parameters*/
	float _minRatio; /*the minimal portion of the query in the optimal local alignment*/
	float _minIdentity; /*the minimal identity in the optimal local alignment*/
	float _gminIdentity; /*minimal identity in the optimal semiglobal alignment*/
	float _gminIdentityRealign;
	uint32_t _maxSeedOcc;
	uint32_t _mapRegionSizeFactor; /*the factor for maximal mapping region*/
	int _minAlignScore;
	int _mapQualReliable;
	int _maxGapSize;
	int64_t _maxDistance;
	int32_t _maxEditDistance; /*maximum edit distance*/
	uint32_t _maxMultiAligns;
	uint32_t _maxSeedPairs; /*maximal number of seed pairs*/
	int _alignType; /*alignment type*/
	int _pairingMode; /*paired-end or mate-paired reads*/
	int _matchScore; /*alignment score for a match*/
	int _seAlign; /*single-end alignment?*/
	float _minPairPriority;

	/*banded global and local aligner*/
	Aligner* _aligner; /*aligner with the user-specified scoring scheme*/
	uint32_t _targetSize;
	uint8_t* _target;

	/*SAM output*/
	SAM* _sam;

	/*processing*/
	vector<uint4> _ranges, _ranges2;
	vector<uint4> _rranges, _rranges2;
	vector<Seed> _seeds, _seeds1, _seeds2;
	vector<Seed> _qualSeeds1, _qualSeeds2;
	vector<int32_t> _bestHits1, _bestHits2;
	vector<SeedPair> _seedPairHeap; /*priority heap*/
	SeedPairComp _seedPairComp;
	set<int64_t> _mapPositions;
	vector<MappingPair> _mappingPairHeap;/*priority heap*/
	MappingPairComp _mappingPairComp;

	/*for SSE2*/
	vector<uint8_t> _sequences;
	vector<int32_t> _seqOffsets;
	vector<AlignScore> _alignScores;

	/*for color-space conversion*/
	vector<uint8_t> _csBuffer;

	/*static variables*/
	static const int COLOR_MM;
	static const int NUCL_MM;
	static const int _bs2csTable[16];
	static const int _cs2bsTable[5][5];
	static const uint8_t _complements[5];

	/*boolean values*/
	bool _colorspace; /*is colorspace*/
	bool _maskAmbiguous; /*mask ambiguous bases in the reference*/
	bool _sensitivePairing;

protected:
};

#endif /* MEMENGINE_H_ */
