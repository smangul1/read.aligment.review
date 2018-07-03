/*
 * MemEngine.cpp
 *
 *  Created on: Dec 24, 2011
 *      Author: yongchao
 */

#include "MemEngine.h"
#include "Bitmap.h"

const int MemEngine::_bs2csTable[] = { 4, 0, 0, 1, 0, 2, 3, 4, 0, 3, 2, 4, 1, 4,
		4, 4 };

/*assume that AN and NA will convert to A for each A*/
const int MemEngine::_cs2bsTable[5][5] = { { 0, 1, 2, 3, 0 }, /*A*/
{ 1, 0, 3, 2, 1 }, /*C*/
{ 2, 3, 1, 0, 2 }, /*G*/
{ 3, 2, 0, 1, 3 }, /*T*/
{ 0, 1, 2, 3, 4 } /*N*/
};
const uint8_t MemEngine::_complements[5] = { 3, 2, 1, 0, 4 };
const int MemEngine::COLOR_MM = 19;
const int MemEngine::NUCL_MM = 25;

MemEngine::MemEngine(Options* options, Genome* rgenome, SAM* sam) {
	_options = options;
	_rgenome = rgenome;
	_sam = sam;
	_rbwt = _rgenome->getBWT();
	_rsa = _rgenome->getSuffixArray();

	/*get the reference sequence length*/
	_bwtSeqLength = _rbwt->getBwtSeqLength();

	/*get the baseBitmap*/
	_baseBitmap = _rgenome->getBaseBitmap();

	/*get the packed genome*/
	_pacGenome = _rgenome->getPacGenome();
	_basePacGenome = _rgenome->getBasePacGenome();

	/*get the minimal seed size*/
	_minIdentity = _options->getMinIdentity() * 100;
	_gminIdentity = _options->getGMinIdentity() * 100;
	_gminIdentityRealign = 0;
	_minRatio = _options->getMinRatio(); /*the minimal portion of the query in the optimal local alignment*/
	_maxSeedOcc = _options->getMaxSeedOcc();
	_minAlignScore = _options->getMinAlignScore();
	_mapQualReliable = 20;
	_mapRegionSizeFactor = 2;
	_maxGapSize = 10;
	_maxEditDistance = _options->getMaxEditDistance();
	_maxMultiAligns = _options->getMaxMultiAligns();
	_maxSeedPairs = _options->getMaxSeedPairs();
	_alignType = _options->getAlignType();
	_matchScore = _options->getMatch();
	_pairingMode = _options->getPairingMode();
	_colorspace = _options->isColorSpace();
	_maskAmbiguous = _options->maskAmbiguous();
	_seAlign = !_options->isPaired();
	_sensitivePairing = _options->isSensitiveParing();
	_minPairPriority = 0.7;

	/*for processing*/
	updateDistance();

	/*create default local aligner*/
	_aligner = new Aligner(options);

	_targetSize = 1024;
	_target = new uint8_t[_targetSize];
	if (_target == NULL) {
		Utils::exit("Memory allocation failed in function %s line %d\n",
				__FUNCTION__, __LINE__);
	}

	/*reserve space*/
	_bestHits1.reserve(128);
	_bestHits2.reserve(128);
}
MemEngine::~MemEngine() {

	/*release the aligner*/
	delete _aligner;
	if (_target) {
		delete[] _target;
	}
}
/*single-end alignment*/
bool MemEngine::align(Sequence& seq, Mapping*& bestMapping, int seedType) {
	/*initialize the mapping*/
	bestMapping = NULL;

	/*get all seeds*/
	_locateSeeds(seq, _seeds, MAXIMAL_EXACT_MATCH_SEED);

	/*get the single-end alignment*/
	if (_seeds.size() > 0) {
		bestMapping = _getAlignment(_aligner, seq, _seeds.size(), &_seeds[0]);
	}

	if (!bestMapping && seedType == MAXIMAL_EXACT_MATCH_SEED) {
		return align(seq, bestMapping, FIXED_LENGTH_STRIPED_SEED);
	}

	return bestMapping != NULL;
}

/*for multiple single-end alignment*/
bool MemEngine::align(Sequence& seq, vector<Mapping*>& mappings, int seedType) {
	pair<set<int64_t>::iterator, bool> ret;
	Mapping* bestMapping;
	Seed* seeds;
	size_t numTopHits;
	int mapQual;

	/*clear the set*/
	_mapPositions.clear();

	/*get all seeds*/
	_locateSeeds(seq, _seeds1, seedType);
	seeds = &_seeds1[0];

	/*get the best hits for the sequence*/
	_getBestHits(_aligner, seq, _seeds1.size(), seeds, _bestHits1, numTopHits,
			mapQual);

	/*get the alignment for the sequence*/
	for (size_t i = 0; i < _bestHits1.size(); ++i) {

		/*get the alignment*/
		bestMapping = _getAlignment(_aligner, seq, seeds[_bestHits1[i]],
				mapQual);

		if (!bestMapping) {
			break;
		}
		/*insert the mapping position to the set to check the uniquness*/
		ret = _mapPositions.insert(bestMapping->_gposition);
		if (ret.second) { /*a new position*/
			mappings.push_back(bestMapping);
		}

		if (mappings.size() >= _maxMultiAligns) {
			break;
		}
	}
	if (mappings.size() == 0 && seedType == MAXIMAL_EXACT_MATCH_SEED) {
		return align(seq, mappings, FIXED_LENGTH_STRIPED_SEED);
	}

	return mappings.size() > 0;
}
/*paired-end alignment, returning multiple equalivalent "best" alignments if applicable*/
bool MemEngine::align(Sequence& seq1, Sequence& seq2,
		vector<Mapping *>& mappings1, vector<Mapping *>& mappings2) {
	int mapQual1, mapQual2;
	size_t numTopHits1, numTopHits2;
	Seed *seeds1, *seeds2;
	int genomeIndex;
	int window, numErrors;
	int64_t genomeStart1 = 0, genomeEnd1 = 0, genomeStart2 = 0, genomeEnd2 = 0;
	int64_t qleft, qright;
	Mapping *bestMapping1, *bestMapping2;
	Mapping *mapping1, *mapping2;
	uint8_t hasRescued = 0;
	int64_t insertSize;
	size_t numSeedPairs;
	pair<set<int64_t>::iterator, bool> ret;
	bool good1, good2, recomp1, recomp2;

	/*initialize*/
	_qualSeeds1.clear();
	_qualSeeds2.clear();
	_mappingPairHeap.clear();
	_mapPositions.clear();
	bestMapping1 = bestMapping2 = NULL;

	/*get seeds for both ends*/
	_locateSeeds(seq1, _seeds1, MAXIMAL_EXACT_MATCH_SEED);
	_locateSeeds(seq2, _seeds2, MAXIMAL_EXACT_MATCH_SEED);

	/*get the best hits for each end*/
	_getBestHits(_aligner, seq1, _seeds1.size(), &_seeds1[0], _bestHits1,
			numTopHits1, mapQual1);
	_getBestHits(_aligner, seq2, _seeds2.size(), &_seeds2[0], _bestHits2,
			numTopHits2, mapQual2);

	/*extract the qualified seeds*/
	_getQualifiedSeeds(_qualSeeds1, _seeds1, _bestHits1, false);
	_getQualifiedSeeds(_qualSeeds2, _seeds2, _bestHits2, false);

	/*get the alignment for the left end*/
	if (_bestHits1.size() > 0) {
		bestMapping1 = _getAlignment(_aligner, seq1, _qualSeeds1[0], mapQual1);
	}

	/*get the alignment for the right end*/
	if (_bestHits2.size() > 0) {
		bestMapping2 = _getAlignment(_aligner, seq2, _qualSeeds2[0], mapQual2);
	}

	if (bestMapping1 == NULL) {
		/*if there is not any qualified hit*/
		_locateSeeds(seq1, _seeds1, FIXED_LENGTH_STRIPED_SEED);

		/*re-calculat the best hits*/
		_getBestHits(_aligner, seq1, _seeds1.size(), &_seeds1[0], _bestHits1,
				numTopHits1, mapQual1);

		/*extract the qualified seeds*/
		_getQualifiedSeeds(_qualSeeds1, _seeds1, _bestHits1, false);

		if (_bestHits1.size() > 0) {
			bestMapping1 = _getAlignment(_aligner, seq1, _qualSeeds1[0],
					mapQual1);
		}
		hasRescued |= 1;
	}

	if (bestMapping2 == NULL) {
		/*if there is not any qualified hit*/
		_locateSeeds(seq2, _seeds2, FIXED_LENGTH_STRIPED_SEED);

		/*re-calculat the best hits*/
		_getBestHits(_aligner, seq2, _seeds2.size(), &_seeds2[0], _bestHits2,
				numTopHits2, mapQual2);

		/*extract the qualified seeds*/
		_getQualifiedSeeds(_qualSeeds2, _seeds2, _bestHits2, false);

		if (_bestHits2.size() > 0) {
			bestMapping2 = _getAlignment(_aligner, seq2, _qualSeeds2[0],
					mapQual2);
		}
		hasRescued |= 2;
	}
	/*get the pointers*/
	seeds1 = &_qualSeeds1[0];
	seeds2 = &_qualSeeds2[0];

	/*check the distance between the sequences*/
	good1 = bestMapping1 && bestMapping1->_mapQual >= _mapQualReliable;
	good2 = bestMapping2 && bestMapping2->_mapQual >= _mapQualReliable;

	if (good1 && good2
			&& (_pairingMode == MATE_PAIRED_READS ?
					bestMapping1->_strand == bestMapping2->_strand :
					bestMapping1->_strand != bestMapping2->_strand)) {
		insertSize = labs(bestMapping1->_position - bestMapping2->_position);
		if (_pairingMode == PAIRED_END_READS) {
			insertSize += bestMapping1->_strand ? seq1._length : seq2._length;
		}
		if (insertSize <= _maxDistance) {
			ret = _mapPositions.insert(bestMapping1->_gposition);
			if (ret.second) {
				/*a new position for the left end*/
				float score = _getHeapPriority(seq1._length, seq2._length,
						bestMapping1->_align->getAlignScore(),
						bestMapping2->_align->getAlignScore());
				_mappingPairHeap.push_back(
						MappingPair(bestMapping1, bestMapping2, score));
				push_heap(_mappingPairHeap.begin(), _mappingPairHeap.end(),
						_mappingPairComp);

				bestMapping1 = bestMapping2 = NULL;
			}
		}
	}

	/*if having found any paired-end alignment*/
	//Utils::log("line: %d #alignes: %ld\n", __LINE__, _mappingPairHeap.size());
	if (_mappingPairHeap.size() >= _maxMultiAligns) {
		/*save the results to the output*/
		size_t numAligns = _mappingPairHeap.size();
		if (numAligns > _maxMultiAligns) {
			numAligns = _maxMultiAligns;
		}
		for (size_t i = 0; i < numAligns; ++i) {
			/*get a new mapping pair*/
			MappingPair& pair = _mappingPairHeap[0];
			mappings1.push_back(pair._left);
			mappings2.push_back(pair._right);

			/*pop the highest-priority element*/
			pop_heap(_mappingPairHeap.begin(), _mappingPairHeap.end(),
					_mappingPairComp);
			_mappingPairHeap.pop_back();
		}
		return true;
	}

	/*attempt to pair reads from seed pairing heuristic*/
	recomp1 = recomp2 = false;
	numErrors = _options->getNumErrors(seq1._length);
	do {
		numSeedPairs = 0;
		_seedPairHeap.clear();
		for (size_t i = 0; i < _qualSeeds1.size(); ++i) {
			Seed& left = seeds1[i];

			/*get the genome index*/
			genomeIndex = left._genomeIndex;

			/*estimate the position of the 5' end*/
			qleft = (left._strand == 0) ?
					left._queryPosition :
					seq1._length - left._queryPosition - 1;
			for (size_t j = 0; j < _qualSeeds2.size(); ++j) {
				Seed& right = seeds2[j];

				/*check the genome index and the sequence strand*/
				if (genomeIndex != right._genomeIndex
						|| (_pairingMode == MATE_PAIRED_READS ?
								left._strand != right._strand :
								left._strand == right._strand)) {
					continue;
				}
				/*estimate the position of the 5'end*/
				qright =
						(right._strand == 0) ?
								right._queryPosition :
								seq2._length - right._queryPosition - 1;

				/*estimate the distance between the two reads (actually, it is the estimated distance of the two 5'ends)*/
				insertSize = labs(
						(int64_t) left._targetPosition - qleft
								- right._targetPosition + qright);
				if (insertSize <= _maxDistance + 2 * numErrors) {
					float score = _getHeapPriority(seq1._length, seq2._length,
							left._alignScore, right._alignScore);
					if (score >= _minPairPriority) {
						_seedPairHeap.push_back(SeedPair(left, right, score));
						push_heap(_seedPairHeap.begin(), _seedPairHeap.end(),
								_seedPairComp);

						if (++numSeedPairs >= _maxSeedPairs) {
							/*reaching the maximum number of seed pairs*/
							break;
						}
					}
				}
			}
			/*reaching the maximum number of seed pairs*/
			if (numSeedPairs >= _maxSeedPairs) {
				break;
			}
		}

		/*pairing reads from paired seed*/
		//Utils::log("#readParis: %ld %ld\n", numSeedPairs, seedPairHeap.size());
		while (_seedPairHeap.size()) {
			/*get a new seed pair*/
			SeedPair seedPair = _seedPairHeap[0];
			/*pop the highest-priority element*/
			pop_heap(_seedPairHeap.begin(), _seedPairHeap.end(), _seedPairComp);
			_seedPairHeap.pop_back();

			/*Utils::log("score: %g pairs: %u %u\n", seedPair._score,
			 seedPair._left._targetPosition, seedPair._right._targetPosition);*/

			/*compute the alignment*/
			mapping1 = _getAlignment(_aligner, seq1, seedPair._left, mapQual1);
			if (!mapping1) {
				continue;
			}
			mapping2 = _getAlignment(_aligner, seq2, seedPair._right, mapQual2);
			if (!mapping2) {
				delete mapping1;
				continue;
			}

			/*calculate the insert size*/
			//check the final mapping positions of the alignment*/
			insertSize = labs(mapping1->_position - mapping2->_position);
			if (_pairingMode == PAIRED_END_READS) {
				insertSize += mapping1->_strand ? seq2._length : seq1._length;
			}
			if (insertSize <= _maxDistance) {
				ret = _mapPositions.insert(mapping1->_gposition);
				if (ret.second) {
					/*a new position for the left end*/
					float score = _getHeapPriority(seq1._length, seq2._length,
							mapping1->_align->getAlignScore(),
							mapping2->_align->getAlignScore());
					_mappingPairHeap.push_back(
							MappingPair(mapping1, mapping2, score));
					push_heap(_mappingPairHeap.begin(), _mappingPairHeap.end(),
							_mappingPairComp);

					mapping1 = mapping2 = NULL;

					if (_mappingPairHeap.size() >= _maxMultiAligns) {
						break;
					}
				}
			}
			if (mapping1)
				delete mapping1;
			if (mapping2)
				delete mapping2;
		}

		/*check if an end has not been rescued*/
		if (hasRescued == 3
				|| _mappingPairHeap.size()
						>= (_sensitivePairing ? _maxMultiAligns : 1)) {
			break;
		}
		if ((hasRescued & 1) == 0) {
			/*if there is not any qualified hit*/
			_locateSeeds(seq1, _seeds1, FIXED_LENGTH_STRIPED_SEED);

			/*re-calculat the best hits*/
			_getBestHits(_aligner, seq1, _seeds1.size(), &_seeds1[0],
					_bestHits1, numTopHits1, mapQual1);

			/*extract the qualified seeds*/
			_getQualifiedSeeds(_qualSeeds1, _seeds1, _bestHits1, true);

			recomp1 = true;
			hasRescued |= 1;
			seeds1 = &_qualSeeds1[0];
		} else if ((hasRescued & 2) == 0) {
			_locateSeeds(seq2, _seeds2, FIXED_LENGTH_STRIPED_SEED);

			/*re-calculat the best hits*/
			_getBestHits(_aligner, seq2, _seeds2.size(), &_seeds2[0],
					_bestHits2, numTopHits2, mapQual2);

			/*extract the qualified seeds*/
			_getQualifiedSeeds(_qualSeeds2, _seeds2, _bestHits2, true);

			recomp2 = true;
			hasRescued |= 2;
			seeds2 = &_qualSeeds2[0];
		}
	} while (1);

	/*if having found any paired-end alignment*/
	//Utils::log("line: %d #alignes: %ld\n", __LINE__, _mappingPairHeap.size());
	if (_mappingPairHeap.size() > 0) {
		/*save the results to the output*/
		size_t numAligns = _mappingPairHeap.size();
		if (numAligns > _maxMultiAligns) {
			numAligns = _maxMultiAligns;
		}
		for (size_t i = 0; i < numAligns; ++i) {
			/*get a new mapping pair*/
			MappingPair& pair = _mappingPairHeap[0];
			mappings1.push_back(pair._left);
			mappings2.push_back(pair._right);

			/*pop the highest-priority element*/
			pop_heap(_mappingPairHeap.begin(), _mappingPairHeap.end(),
					_mappingPairComp);
			_mappingPairHeap.pop_back();
		}
		return true;
	}

	/*recompute the best alignments*/
	if (recomp1 && _bestHits1.size() > 0) {
		if (bestMapping1) {
			delete bestMapping1;
		}
		bestMapping1 = _getAlignment(_aligner, seq1, _qualSeeds1[0], mapQual1);
	}

	if (recomp2 && _bestHits2.size() > 0) {
		if (bestMapping2) {
			delete bestMapping2;
		}
		bestMapping2 = _getAlignment(_aligner, seq2, _qualSeeds2[0], mapQual2);
	}

	if (recomp1 || recomp2) {
		/*check the distance between the sequences*/
		good1 = bestMapping1 && bestMapping1->_mapQual >= _mapQualReliable;
		good2 = bestMapping2 && bestMapping2->_mapQual >= _mapQualReliable;
		if (good1 && good2
				&& (_pairingMode == MATE_PAIRED_READS ?
						bestMapping1->_strand == bestMapping2->_strand :
						bestMapping1->_strand != bestMapping2->_strand)) {
			insertSize = labs(
					bestMapping1->_position - bestMapping2->_position);
			if (_pairingMode == PAIRED_END_READS) {
				insertSize +=
						bestMapping1->_strand ? seq1._length : seq2._length;
			}
			if (insertSize <= _maxDistance) {
				//Utils::log("line: %d \n", __LINE__);
				/*the two reads are paired and return*/
				mappings1.push_back(bestMapping1);
				mappings2.push_back(bestMapping2);

				return true;
			}
		}
	}

	/*try to find the best alignment for the right sequence*/
	Mapping *mateMapping1 = NULL, *mateMapping2 = NULL;
	if (_options->rescueMate()) {
		if (good1) {
			//Utils::log("map1: %u %u\n", bestMapping1->_position, bestMapping1->_strand);
			window = _options->getNumErrors(seq2._length) * 2;

			/*calcualte the mapping region of the mate read*/
			_getMateRegion(genomeStart1, genomeEnd1, bestMapping1->_gposition,
					bestMapping1->_strand, seq1._length, seq2._length);

			/*refine the region*/
			_rgenome->refineRegionRange(bestMapping1->_genomeIndex,
					genomeStart1, genomeEnd1);

			/*perform the alignment*/
			mateMapping1 = _getAlignment(_aligner, seq2,
					_pairingMode == MATE_PAIRED_READS ?
							bestMapping1->_strand : 1 - bestMapping1->_strand,
					window, bestMapping1->_genomeIndex, genomeStart1,
					genomeEnd1 - genomeStart1 + 1, mapQual1, _colorspace ? true : false);

			/*succeeded in finding an alignment*/
			if (mateMapping1) {
				mappings1.push_back(bestMapping1);
				mappings2.push_back(mateMapping1);
				bestMapping1 = NULL;

				if (bestMapping2) {
					delete bestMapping2;
				}
				return true;
			}
		}

		/*rescue the alignment from the right sequence*/
		if (good2) {
			window = _options->getNumErrors(seq1._length) * 2;
			/*calcualte the mapping region of the mate read*/
			_getMateRegion(genomeStart2, genomeEnd2, bestMapping2->_gposition,
					bestMapping2->_strand, seq2._length, seq1._length);

			/*refine the region*/
			_rgenome->refineRegionRange(bestMapping2->_genomeIndex,
					genomeStart2, genomeEnd2);

			/*perform the alignment*/
			mateMapping2 = _getAlignment(_aligner, seq1,
					_pairingMode == MATE_PAIRED_READS ?
							bestMapping2->_strand : 1 - bestMapping2->_strand,
					window, bestMapping2->_genomeIndex, genomeStart2,
					genomeEnd2 - genomeStart2 + 1, mapQual2, _colorspace ? true : false);

			/*failed to find an alignment*/
			if (mateMapping2) {
				/*output the paired-end alignments*/
				mappings1.push_back(mateMapping2);
				mappings2.push_back(bestMapping2);
				bestMapping2 = NULL;

				if (bestMapping1) {
					delete bestMapping1;
				}
				return true;
			}
		}
	}

	/*twice rescuing*/
	size_t numRescues = _options->rescueTwice();
	if (numRescues > 0) {
		if (!good1 && bestMapping1) {
			window = _options->getNumErrors(seq2._length) * 2;
			/*perform the alignment*/
			mateMapping1 = _getAlignment(_aligner, seq2,
					_pairingMode == MATE_PAIRED_READS ?
							bestMapping1->_strand : 1 - bestMapping1->_strand,
					window, bestMapping1->_genomeIndex, genomeStart1,
					genomeEnd1 - genomeStart1 + 1,
					SW_MAP_QUALITY_SCORE * mapQual1, true);

			/*succeeded in finding an alignment*/
			if (mateMapping1) {
				/*output the paired-end alignments*/
				mappings1.push_back(bestMapping1);
				mappings2.push_back(mateMapping1);
				bestMapping1 = NULL;
				if (bestMapping2) {
					delete bestMapping2;
				}
				return true;
			}
		}
		if (!good2 && bestMapping2) {
			window = _options->getNumErrors(seq1._length) * 2;
			/*perform the alignment*/
			mateMapping2 = _getAlignment(_aligner, seq1,
					_pairingMode == MATE_PAIRED_READS ?
							bestMapping2->_strand : 1 - bestMapping2->_strand,
					window, bestMapping2->_genomeIndex, genomeStart2,
					genomeEnd2 - genomeStart2 + 1,
					SW_MAP_QUALITY_SCORE * mapQual2, true);

			/*failed to find an alignment*/
			if (mateMapping2) {
				/*output the paired-end alignments*/
				mappings1.push_back(mateMapping2);
				mappings2.push_back(bestMapping2);
				bestMapping2 = NULL;
				if (bestMapping1) {
					delete bestMapping1;
				}
				return true;
			}
		}
		/*rescure from other top hits*/
		mapping1 = NULL;
		numTopHits1 = min(numRescues, _qualSeeds1.size());
		for (size_t index = 1; index < numTopHits1; ++index) {
			Seed& seed = seeds1[index];

			/*generate an alignment*/
			if ((mapping1 = _getAlignment(_aligner, seq1, seed, mapQual1))
					== NULL) {
				break;
			}

			/*rescue from the current alignment*/
			window = _options->getNumErrors(seq2._length) * 2;
			/*perform the alignment*/
			//Utils::log("strand: %d genomeStart1: %ld genomeEnd: %ld\n",  _pairingMode == MATE_PAIRED_READS ? bestMapping1->_strand : 1 - bestMapping1->_strand, genomeStart1, genomeEnd1);
			/*calcualte the mapping region of the mate read*/
			_getMateRegion(genomeStart1, genomeEnd1, mapping1->_gposition,
					mapping1->_strand, seq1._length, seq2._length);

			/*refine the region*/
			_rgenome->refineRegionRange(mapping1->_genomeIndex, genomeStart1,
					genomeEnd1);

			mateMapping1 = _getAlignment(_aligner, seq2,
					_pairingMode == MATE_PAIRED_READS ?
							mapping1->_strand : 1 - mapping1->_strand, window,
					mapping1->_genomeIndex, genomeStart1,
					genomeEnd1 - genomeStart1 + 1,
					SW_MAP_QUALITY_SCORE * mapQual1, true);

			/*succeeded in finding an alignment*/
			if (mateMapping1) {
				/*output the paired-end alignments*/
				mappings1.push_back(mapping1);
				mappings2.push_back(mateMapping1);
				if (bestMapping1) {
					delete bestMapping1;
				}
				if (bestMapping2) {
					delete bestMapping2;
				}
				return true;
			} else {
				delete mapping1;
			}
		}

		/*find the next best with different strand*/
		mapping2 = NULL;
		numTopHits2 = min(numTopHits2, _qualSeeds2.size());
		for (size_t index = 1; index < numTopHits2; ++index) {
			Seed& seed = seeds2[index];

			/*get a new alignment*/
			if ((mapping2 = _getAlignment(_aligner, seq2, seed, mapQual2))
					== NULL) {
				break;
			}

			/*rescue the alignment*/
			window = _options->getNumErrors(seq1._length) * 2;
			/*calcualte the mapping region of the mate read*/
			_getMateRegion(genomeStart2, genomeEnd2, mapping2->_gposition,
					mapping2->_strand, seq2._length, seq1._length);

			/*refine the region*/
			_rgenome->refineRegionRange(mapping2->_genomeIndex, genomeStart2,
					genomeEnd2);

			/*perform the alignment*/
			mateMapping2 = _getAlignment(_aligner, seq1,
					_pairingMode == MATE_PAIRED_READS ?
							mapping2->_strand : 1 - mapping2->_strand, window,
					mapping2->_genomeIndex, genomeStart2,
					genomeEnd2 - genomeStart2 + 1,
					SW_MAP_QUALITY_SCORE * mapQual2, true);

			/*failed to find an alignment*/
			if (mateMapping2) {
				/*output the paired-end alignments*/
				mappings1.push_back(mateMapping2);
				mappings2.push_back(mapping2);

				if (bestMapping1) {
					delete bestMapping1;
				}
				if (bestMapping2) {
					delete bestMapping2;
				}
				return true;
			} else {
				delete mapping2;
			}
		}
	}

	/*If failed to find any paired-end alignment, will output the best alignments as single-end ones*/
	if (bestMapping1) {
		mappings1.push_back(bestMapping1);
	}
	if (bestMapping2) {
		mappings2.push_back(bestMapping2);
	}
	return false;
}

int32_t MemEngine::_getBestHit(Aligner* aligner, Sequence& seq, size_t numSeeds,
		Seed* seeds, int& mapQual) {
	int32_t bestScore, bestScore2, bestSeedIndex;
	size_t totalNseeds, seedIndex;
	bool found = false;
	AlignScore* alignScores;
	int32_t numErrors = _options->getNumErrors(seq._length);
	int64_t targetPosition, targetPosition2;

	/*get the alignment scores*/
	totalNseeds = _getAlignmentScores(aligner, seq, numSeeds, seeds);

	mapQual = 0;
//Utils::log("totalNSeeds: %d\n", totalNseeds);
	if (totalNseeds > 0) {
		/*select the best alignments*/
		stable_sort(_alignScores.begin(), _alignScores.end());
		alignScores = &_alignScores[0];
#if 0
		/*print out the alignment scores*/
		Utils::log("-------------Ranked seeds:----------------\n");
		for (size_t i = 0; i < totalNseeds; ++i) {
			Seed* p = &seeds[alignScores[i]._seedIndex];
			fprintf(stderr, "%u %u %u %u %u\n", p->_queryPosition,
					p->_targetPosition, p->_seedLength, p->_strand,
					alignScores[i]._score);
		}
#endif
		bestScore = alignScores->_score;
		if (bestScore < _minAlignScore) {
			return -1;
		}
		bestSeedIndex = alignScores->_seedIndex;

		//calculate mapping quality score
		//estimate the starting position of the alignment
		targetPosition = alignScores->_targetPosition;
		targetPosition -= alignScores->_seedPosition;
		++alignScores;
		for (seedIndex = 1; seedIndex < totalNseeds; ++seedIndex) {
			/*get the next best score*/
			bestScore2 = alignScores->_score;
			if (bestScore2 != bestScore) {
				found = true;
				break;
			}

			/*estimate the starting position of the alignment*/
			targetPosition2 = alignScores->_targetPosition;
			targetPosition2 -= alignScores->_seedPosition;
			//Utils::log("pos1 %ld pos2 %ld\n", targetPosition, targetPosition2);
			if (bestScore2 == bestScore
					&& labs(targetPosition - targetPosition2)
							>= _maxGapSize + numErrors) {
				found = true;
				break;
			}
			++alignScores;
		}
		mapQual =
				(found == true) ?
						DEFAULT_MAX_MAP_QUAL * (bestScore - bestScore2)
								/ bestScore :
						DEFAULT_MAX_MAP_QUAL;

#if 0
		/*print out the alignment scores*/
		Utils::log("-------------Seed Score:----------------\n");
		Utils::log("mapping quality: %d\n", mapQual);
		alignScores = &_alignScores[0];
		for (size_t i = 0; i < totalNseeds; ++i) {
			Utils::log("query %d strand %d score %d target %u\n",
					seeds[alignScores[i]._seedIndex]._queryPosition,
					seeds[alignScores[i]._seedIndex]._strand,
					alignScores[i]._score, alignScores[i]._targetPosition);
		}
#endif

		return bestSeedIndex;
	}
	return -1;
}

/*for paired-end alignment*/
int32_t MemEngine::_getBestHits(Aligner* aligner, Sequence& seq,
		size_t numSeeds, Seed* seeds, vector<int32_t>& bestHits,
		size_t& numTopSeeds, int& mapQual) {
	size_t totalNseeds, seedIndex;
	int32_t bestScore, bestScore2, bestScoreDiff;
	bool found = false;
	int64_t targetPosition, targetPosition2;
	AlignScore* alignScores;
	int32_t numErrors = _options->getNumErrors(seq._length);

	/*get the alignment scores*/
	totalNseeds = _getAlignmentScores(aligner, seq, numSeeds, seeds);

	mapQual = 0;
	numTopSeeds = 0;
	bestHits.clear();

	/*calculate the minimal alignment score for a qualified seed*/
	if (totalNseeds > 0) {
		/*select the best alignments*/
		stable_sort(_alignScores.begin(), _alignScores.end());
		alignScores = &_alignScores[0];
#if 0
		/*print out the alignment scores*/
		Utils::log("-------------Ranked seeds:----------------\n");
		seq.print(stderr);
		for (size_t i = 0; i < totalNseeds; ++i) {
			Seed* p = &seeds[alignScores[i]._seedIndex];
			fprintf(stderr, "%u %u %u %u %u\n", p->_queryPosition,
					p->_targetPosition, p->_seedLength, p->_strand,
					alignScores[i]._score);
		}
#endif

		/*select the best hits*/
		bestScore = alignScores->_score;

		/*check the minimal alignment score*/
		if (bestScore < _minAlignScore) {
			return 0;
		}
		//select the top seeds
		bestScoreDiff = bestScore;
		//estimate the starting position of the alignment
		targetPosition = alignScores->_targetPosition;
		targetPosition -= alignScores->_seedPosition;
		for (seedIndex = 0; seedIndex < totalNseeds; ++seedIndex) {

			/*get the next best score*/
			bestScore2 = alignScores->_score;
			if (bestScore2 != bestScore) {
				if (!found) {
					found = true;
					/*optimal local alignment score diff*/
					bestScoreDiff = bestScore - bestScore2;
				}
				break;
			}
			/*estimate the starting position of the alignment*/
			targetPosition2 = alignScores->_targetPosition;
			targetPosition2 -= alignScores->_seedPosition;
			if (!found
					&& (bestScore2 == bestScore
							&& labs(targetPosition - targetPosition2)
									>= _maxGapSize + numErrors)) {
				found = true;
				/*optimal local aignment score diff*/
				bestScoreDiff = 0;
			}
			bestHits.push_back(alignScores->_seedIndex);

			/*save the alignment score*/
			seeds[alignScores->_seedIndex]._alignScore = alignScores->_score;
			++alignScores;
		}
		//calculate mapping quality score
		mapQual =
				(found == true) ?
						DEFAULT_MAX_MAP_QUAL * bestScoreDiff / bestScore :
						DEFAULT_MAX_MAP_QUAL;

		/*select second-best top hits*/
		numTopSeeds = bestHits.size();

		/*fill the best hits*/
		alignScores = &_alignScores[seedIndex];
		for (; seedIndex < totalNseeds; ++seedIndex) {
			if (alignScores->_score < _minAlignScore) {
				break;
			}
			bestHits.push_back(alignScores->_seedIndex);
			/*save the alignment score*/
			seeds[alignScores->_seedIndex]._alignScore = alignScores->_score;
			++alignScores;
		}
	}
#if 0
	/*print out the alignment scores*/
	Utils::log("-------------Top seeds:----------------\n");
	seq.print(stderr);
	alignScores = &_alignScores[0];
	for (size_t i = 0; i < bestHits.size(); ++i) {
		Seed* p = &seeds[bestHits[i]];
		fprintf(stderr, "%u %u %u %u %u (genome)%d\n", p->_queryPosition,
				p->_targetPosition, p->_seedLength, p->_strand,
				alignScores[i]._score, p->_genomeIndex);
	}
#endif
	return numTopSeeds;
}
size_t MemEngine::_getAlignmentScores(Aligner* aligner, Sequence& seq,
		size_t numSeeds, Seed* seeds) {
	int64_t lowerBound, upperBound, genomeLength;
	int32_t nseeds, offset, nseedsPerBatch;
	size_t seedIndex;
	AlignScore* scorePtr, *globalScorePtr;
	uint8_t* sequences;
	int32_t* seqOffsets;
	uint8_t* ptr = NULL;
	int32_t batchSize = _options->getMaxSeedsPerBatch();
	uint8_t* query;

//Utils::log("batchSize: %d\n", batchSize);
	/*if there are no seeds*/
	if (numSeeds == 0) {
		_alignScores.clear();
		return 0;
	}

//Utils::log("numSeeds: %d\n", numSeeds);
	/*alignment score buffer*/
	_alignScores.resize(numSeeds);

	/*resize the sequence buffer*/
	_seqOffsets.resize(batchSize);
	_sequences.resize(
			min((size_t) batchSize, numSeeds)
					* (_mapRegionSizeFactor * (((seq._length + 3) >> 2) << 2)
							+ 1));

	/*calculate alignment score for all seed extensions*/
	nseeds = 0;
	globalScorePtr = &_alignScores[0];
	sequences = &_sequences[0];
	seqOffsets = &_seqOffsets[0];
	/*evaluate all seeds*/
	for (int strand = 0; strand < 2; ++strand) {

		/*get the query sequence*/
		query = (strand == 0) ? seq._bases : seq._rbases;
		for (uint32_t i = 0; i < seq._length; ++i) {
			query[i] += 1;
		}

		/*initialize variables*/
		offset = 0;
		nseedsPerBatch = 0;
		scorePtr = globalScorePtr + nseeds;
		for (seedIndex = 0; seedIndex < numSeeds; ++seedIndex) {
			Seed& seed = seeds[seedIndex];
			//Utils::log("strand: %d, tlength: %d length %d\n", seed._strand, seq._tlength, seq._length);
			if (seed._strand != strand) {
				continue;
			}

			/*perform banded local alignment using the highest-scoring seeds*/
			lowerBound = seed._targetPosition;
			lowerBound -= _mapRegionSizeFactor * (seed._queryPosition + 1);
			upperBound = seed._targetPosition;
			upperBound += seed._seedLength
					+ _mapRegionSizeFactor
							* (seq._length - seed._queryPosition
									- seed._seedLength);

			//Utils::log("[before] lowerBound: %ld upperBound: %ld seedLength %ld\n", lowerBound, upperBound, seed._seedLength);
			/*refine the genome region*/
			/*get the genome index*/
			_rgenome->refineRegionRange(seed._genomeIndex, lowerBound,
					upperBound);

			genomeLength = upperBound - lowerBound + 1;
			if (genomeLength <= 0) {
				continue;
			}
			/*save the offset*/
			seqOffsets[nseedsPerBatch++] = offset;

			/*save the seed index in the global alignment score buffer*/
			globalScorePtr[nseeds]._seedLength = seed._seedLength;
			globalScorePtr[nseeds]._seedPosition = seed._queryPosition;
			globalScorePtr[nseeds]._seedIndex = seedIndex;
			globalScorePtr[nseeds]._targetPosition = seed._targetPosition;
			++nseeds;

			/*load the target sequence slice*/
			/*check the sequence length*/
			ptr = sequences + offset;
			if (!_maskAmbiguous) {
				for (int64_t i = 0, j = lowerBound; i < genomeLength;
						++i, ++j) {
					ptr[i] = ((_pacGenome[j >> 2] >> ((~j & 3) << 1)) & 3) + 1;
				}
			} else {
				for (int64_t i = 0, j = lowerBound; i < genomeLength;
						++i, ++j) {
					ptr[i] =
							((_baseBitmap[j >> 3] >> (j & 7)) & 1) ?
									UNKNOWN_BASE :
									((_pacGenome[j >> 2] >> ((~j & 3) << 1)) & 3)
											+ 1;
				}
			}
			int64_t pGenomeLength = ((genomeLength + 3) >> 2) << 2;
			for (int64_t i = genomeLength; i < pGenomeLength; ++i) {
				ptr[i] = 6;
			}
			ptr[pGenomeLength] = 0; /*separate the sequences*/

			/*increase the offset*/
			offset += pGenomeLength + 1;

			/*check the batch size*/
			if (nseedsPerBatch >= batchSize) {
				//Utils::log("nseedsPerBatch: %d nseeds %d offset: %d\n", nseedsPerBatch, nseeds, offset);
				/*calculate the optimal local alignment scores*/
				aligner->lalignScore(query, seq._length, sequences, seqOffsets,
						nseedsPerBatch, scorePtr);

				/*update other variables*/
				scorePtr += nseedsPerBatch;
				nseedsPerBatch = 0;
				offset = 0;
			}
		}
		/*check the batch size*/
		if (nseedsPerBatch > 0) {
			//Utils::log("nseedsPerBatch: %d nseeds %d offset: %d\n", nseedsPerBatch, nseeds, offset);
			/*calculate the optimal local alignment scores*/
			aligner->lalignScore(query, seq._length, sequences, seqOffsets,
					nseedsPerBatch, scorePtr);
		}

		/*restore the query*/
		for (uint32_t i = 0; i < seq._length; ++i) {
			query[i] -= 1;
		}
	}
	/*set the total number of seeds*/
//Utils::log("nseeds %d numSeeds: %d\n", nseeds, numSeeds);
	_alignScores.resize(nseeds);
	return nseeds;
}
Mapping* MemEngine::_getAlignment(Aligner* aligner, Sequence& seq, Seed& seed,
		int mapQual) {
	int32_t window, genomeIndex;
	vector<CigarAlign*> aligns;
	int64_t lowerBound, upperBound;
	Mapping* mapping = NULL;
	int strand = seed._strand;

	/*get the genome index*/
	genomeIndex = seed._genomeIndex;

	/*check if it is an exact-match*/
	if (seed._seedLength == seq._length) {
		mapping = new Mapping(new CigarAlign(seq, _matchScore),
				seed._targetPosition - _rgenome->getGenomeOffset(genomeIndex)
						+ 1, seed._targetPosition, strand, genomeIndex, mapQual,
				mapQual);
	} else {
		/*calculate the banded window width for banded local/global alignment*/
		window = _options->getNumErrors(seq._length) * 2;

		/*caluculate the mapping region on the genome*/
		upperBound = lowerBound = seed._targetPosition;
		lowerBound -= _mapRegionSizeFactor * (seed._queryPosition + 1);
		upperBound +=
				seed._seedLength
						+ _mapRegionSizeFactor
								* (seq._length - seed._queryPosition
										- seed._seedLength);

		//Utils::log("queryPosition: %d strand: %d targetPosition:%u\n", seed._queryPosition, strand, seed._targetPosition);
		//Utils::log("[before] lowerBound: %ld upperBound: %ld seedLength %ld\n", lowerBound, upperBound, seed._seedLength);
		/*refine the genome region*/
		_rgenome->refineRegionRange(genomeIndex, lowerBound, upperBound);
		//Utils::log("[after] lowerBound: %ld upperBound: %ld seedLength %ld\n", lowerBound, upperBound, seed._seedLength);

		int64_t genomeLength = upperBound - lowerBound + 1;
		if (genomeLength <= 0) {
			return mapping;
		}
		/*get the short read alignment*/
		mapping = _getAlignment(aligner, seq, strand, window, genomeIndex,
				lowerBound, genomeLength, mapQual, false);
	}
	/*convert color-space to base-space*/
	if (_colorspace && mapping) {
		uint8_t* bases;
		Mapping* mapping2;
		/*convert read first*/
		pair<uint8_t*, int32_t> newRead = _cs2bsReadConvertBWA(seq, mapping);

		/*re-alignment using semi-global alignment, if there are indels;
		 * local alignment, otherwise*/
		bases = strand ? newRead.first + newRead.second : newRead.first;
		mapping2 = _getAlignment(aligner, mapping->_align->getAlignType(),
				_basePacGenome, bases, newRead.second, strand,
				mapping->_genomeIndex, mapping->_gposition, mapping->_align,
				mapping->_mapQualBase, _gminIdentityRealign);

		/*reset the mapping*/
		delete mapping;
		mapping = mapping2;

		if (mapping) {
			mapping->setData(newRead.first, newRead.second);
		} else {
			if (newRead.first) {
				delete[] newRead.first;
			}
		}
	}
	/*check the edit distance and minimal mapping quality score*/
	if (mapping && mapping->_align->getEditDistance() > _maxEditDistance) {
		/*release the mapping*/
		delete mapping;
		mapping = NULL;
	}
	return mapping;
}

/*normal single-end alignment*/
Mapping* MemEngine::_getAlignment(Aligner* aligner, Sequence& seq, int strand,
		int window, int genomeIndex, uint32_t genomeStart, size_t genomeLength,
		int mapQual, bool onlyLocalAlign) {
	uint8_t* bases;
	int low, up;
	int diff;
	vector<CigarAlign*> aligns;
	Mapping* mapping = NULL;

	/*check the sequence length*/
	if (seq._tlength < GLOBAL_MIN_SEED_SIZE) {
		return NULL;
	}

	/*check the number of leading and trailing bases*/
	bases = (strand == 0) ? seq._bases : seq._rbases;

//Utils::log("genome start: %u genome length %u\n", genomeStart, genomeLength);
	/*get the target sequence*/
	if (genomeLength >= _targetSize) {
		/*resize the target*/
		_targetSize = genomeLength << 1;
		if (_target) {
			delete[] _target;
		}
		_target = new uint8_t[_targetSize];
		if (_target == NULL) {
			Utils::exit("Memory allocation failed in function %s line %d\n",
					__FUNCTION__, __LINE__);
		}
	}
	/*load the target sequence slice*/
	if (!_maskAmbiguous) {
		for (uint32_t i = 0, j = genomeStart; i < genomeLength; ++i, ++j) {
			_target[i] = (_pacGenome[j >> 2] >> ((~j & 3) << 1)) & 3;
		}
	} else {
		for (uint32_t i = 0, j = genomeStart; i < genomeLength; ++i, ++j) {
			_target[i] =
					((_baseBitmap[j >> 3] >> (j & 7)) & 1) ?
							UNKNOWN_BASE :
							((_pacGenome[j >> 2] >> ((~j & 3) << 1)) & 3);
		}
	}
	/***********************************
	 * perform local alignment
	 * **************************************/
	bool swapped = false;
	if (genomeLength >= seq._length) {
		diff = genomeLength - seq._length;
	} else {
		diff = seq._length - genomeLength;
		swapped = true;
	}
	/*make sure that low <= min(0, -diff) and up >= max(0, -diff)*/
	low = -diff - window / 2;
	up = window + low;
	if (up < window) {
		up = window;
	}
	/*Utils::log("diff: %d window %d low %d up:%d\n", diff, window, low, up);*/

	/*perform banded local alignment and output the alignment for the query*/
	if (swapped == false) {
		/*the second  sequence is the query*/
		aligns = aligner->lalignPath(_target, bases, genomeLength, seq._length,
				low, up, 2);
	} else {

		/*the first sequence is the query*/
		aligns = aligner->lalignPath(bases, _target, seq._length, genomeLength,
				low, up, 1);
	}
	/*failed to find an alignment*/
	if (aligns.size() == 0) {
		return NULL;
	}
	/*Utils::log("start1: %d end1: %d start2: %d end2:%d\n", aligns[0]->getStart(), aligns[0]->getEnd(), aligns[0]->getMateStart(),
	 aligns[0]->getMateEnd());*/

	/*optimal local alignment is already good enough?*/
	if (aligns[0]->getNumBases1() >= _minRatio * seq._length
			&& aligns[0]->getIdentity() >= _minIdentity) {

		/*extend the alignment*/
		aligns[0]->extendCigar(seq._length);

		/*mapping positions on the target sequence, indexed from 1*/
		int64_t gPosition = genomeStart + aligns[0]->getMateStart();
		int64_t mapPosition = gPosition - _rgenome->getGenomeOffset(genomeIndex)
				+ 1;

		/*it is possible that the mapping position is less than 0 after extending the alignment to the begining*/
		if (mapPosition >= 0) {
			mapping = new Mapping(aligns[0], mapPosition, gPosition, strand,
					genomeIndex,
					mapQual * aligns[0]->getNumBases1() / seq._length, mapQual);
			aligns[0] = NULL;
		} else {
			Utils::log("negative mapping position\n");
		}
	} else if (!onlyLocalAlign && _alignType == HAVE_LOCAL_SEMI_GLOBAL_ALIGN) {
		/*perform semiglobal alignment around the optimal local alignment*/
		mapping = _getAlignment(aligner, REALIGN_TYPE_SEMIGLOBAL, _pacGenome,
				strand ? seq._rbases : seq._bases, seq._length, strand,
				genomeIndex, genomeStart + aligns[0]->getMateStart(), aligns[0],
				mapQual, _gminIdentity);
	}
	/*release the old alignments*/
	for (size_t i = 0; i < aligns.size(); ++i) {
		if (aligns[i]) {
			delete aligns[i];
		}
	}
	aligns.clear();

	/*paired-end alignment read rescue for color-space reads*/
	if (onlyLocalAlign) {
		uint8_t* bases;
		Mapping* mapping2;
		if (_colorspace && mapping) {
			pair<uint8_t*, int32_t> newRead = _cs2bsReadConvertBWA(seq,
					mapping);

			/*re-alignment using local alignment*/
			bases = strand ? newRead.first + newRead.second : newRead.first;
			mapping2 = _getAlignment(aligner, REALIGN_TYPE_LOCAL,
					_basePacGenome, bases, newRead.second, strand,
					mapping->_genomeIndex, mapping->_gposition, mapping->_align,
					mapping->_mapQualBase, _gminIdentityRealign);

			/*reset the mapping*/
			delete mapping;
			mapping = mapping2;

			if (mapping) {
				mapping->setData(newRead.first, newRead.second);
			} else {
				if (newRead.first) {
					delete[] newRead.first;
				}
			}
		}

		/*check the edit distance*/
		if (mapping && mapping->_align->getEditDistance() > _maxEditDistance) {
			/*release the mapping*/
			delete mapping;
			mapping = NULL;
		}
	}

	return mapping;
}

/*this function is mainly for re-alignment using optimal local alignment as seed*/
Mapping* MemEngine::_getAlignment(Aligner* aligner, int alignType,
		uint8_t* genome, uint8_t* seqBases, uint32_t seqLength, int strand,
		int64_t genomeIndex, int64_t targetPosition, CigarAlign* align,
		int mapQual, float minIdentity) {
	/*consider the optimal local alignment as the seed*/
	int64_t lowerBound, upperBound;
	uint32_t genomeLength, genomeStart;
	int64_t queryPosition = align->getStart(); /*the starting position of the local alignment on the query*/
	int64_t seedLength = align->getAlignLength(); /*the length of the optimal local alignment*/
//int64_t nerrors = _options->getNumErrors(seqLength);
	vector<CigarAlign*> aligns;
	Mapping* mapping = NULL;

	/*check the sequence length*/
	if (seqLength < GLOBAL_MIN_SEED_SIZE) {
		return NULL;
	}

	/*caluclate the mapping region around the optimal local alignment*/
	if (alignType == REALIGN_TYPE_LOCAL) {
		upperBound = lowerBound = targetPosition;
		lowerBound -= _mapRegionSizeFactor * queryPosition;
		upperBound += seedLength
				+ _mapRegionSizeFactor
						* max(0,
								(int32_t) (seqLength - queryPosition
										- seedLength)) - 1;
	} else {
		upperBound = lowerBound = targetPosition;
		lowerBound -= queryPosition;
		upperBound += seedLength
				+ max(0, (int32_t) (seqLength - queryPosition - seedLength))
				- 1;
	}

	/*refine the mapping region*/
	_rgenome->refineRegionRange(genomeIndex, lowerBound, upperBound);
	genomeStart = lowerBound;
	genomeLength = upperBound - lowerBound + 1;
	/*Utils::log("genomeStart: %ld genomeLength %ld queryPosition %ld\n",
	 genomeStart, genomeLength, queryPosition);*/

	/*get the sequence fragment*/
	if (genomeLength >= _targetSize) {
		/*resize the target*/
		_targetSize = genomeLength << 1;
		if (_target) {
			delete[] _target;
		}
		_target = new uint8_t[_targetSize];
		if (_target == NULL) {
			Utils::exit("Memory allocation failed in function %s line %d\n",
					__FUNCTION__, __LINE__);
		}
	}
	if (!_maskAmbiguous) {
		for (uint32_t i = 0, j = genomeStart; i < genomeLength; ++i, ++j) {
			_target[i] = (genome[j >> 2] >> ((~j & 3) << 1)) & 3;
		}
	} else {
		for (uint32_t i = 0, j = genomeStart; i < genomeLength; ++i, ++j) {
			_target[i] =
					((_baseBitmap[j >> 3] >> (j & 7)) & 1) ?
							UNKNOWN_BASE :
							((genome[j >> 2] >> ((~j & 3) << 1)) & 3);
		}
	}

	/*perform local/semiglobal alignment*/
	if (alignType == REALIGN_TYPE_LOCAL) {
		/*perform local alignment*/
		bool swapped = false;
		int32_t diff, low, up;
		int32_t window = _options->getNumErrors(seqLength) * 2;
		if (genomeLength >= seqLength) {
			diff = genomeLength - seqLength;
		} else {
			diff = seqLength - genomeLength;
			swapped = true;
		}
		/*make sure that low <= min(0, -diff) and up >= max(0, -diff)*/
		low = -diff - window / 2;
		up = window + low;
		if (up < window) {
			up = window;
		}
		/*Utils::log("diff: %d window %d low %d up:%d\n", diff, window, low, up);*/

		/*perform banded local alignment and output the alignment for the query*/
		if (swapped == false) {
			/*the second  sequence is the query*/
			aligns = aligner->lalignPath(_target, seqBases, genomeLength,
					seqLength, low, up, 2);
		} else {

			/*the first sequence is the query*/
			aligns = aligner->lalignPath(seqBases, _target, seqLength,
					genomeLength, low, up, 1);
		}
	} else {
		aligns = aligner->galignPath(seqBases, _target, seqLength,
				genomeLength);
	}

	/*check the alignment*/
	if (aligns.size() == 0) {
		return NULL;
	}

	/*mapping positions on the target sequence, indexed from 1*/
//Utils::log("identity: %g\n", aligns[0]->getIdentity());
	if (aligns[0]->getNumBases1() >= _minRatio * seqLength
			&& aligns[0]->getIdentity() >= minIdentity) {

		/*extend the alignment*/
		aligns[0]->extendCigar(seqLength);

		/*calculate the mapping positions*/
		int64_t gPosition = genomeStart + aligns[0]->getMateStart();
		int64_t mapPosition = gPosition - _rgenome->getGenomeOffset(genomeIndex)
				+ 1;

		/*it is possible that the mapping position is less than 0 after extending the alignment to the begining*/
		if (mapPosition >= 0) {
			mapping = new Mapping(aligns[0], mapPosition, gPosition, strand,
					genomeIndex,
					mapQual * aligns[0]->getNumBases1() / seqLength, mapQual);
			aligns[0] = NULL;
		} else {
			Utils::log("negative mapping position\n");
		}
	}
	/*release the old alignments*/
	for (size_t i = 0; i < aligns.size(); ++i) {
		if (aligns[i]) {
			delete aligns[i];
		}
	}
	aligns.clear();

	return mapping;
}
uint32_t MemEngine::_genMEMSeeds(Sequence& seq, int strand, uint4* ranges,
		uint32_t minSeedSize) {
	uint8_t ch;
	bool lastStopAtN;
	uint2 range, range2;
	/*For each position of the sequence, we search the maximal exact match*/
	uint32_t startPos, endPos, lastPos, numRanges;
	uint8_t* bases = (strand == 0) ? seq._bases : seq._rbases;
	uint32_t numSuffixes = seq._length - minSeedSize + 1;

	startPos = 0;
	lastPos = 0;
	numRanges = 0;
	while (startPos < numSuffixes) {
		//for each starting position in the query
		lastStopAtN = false;
		range = make_uint2(0, _bwtSeqLength);
		range2 = make_uint2(1, 0);
		for (endPos = startPos; endPos < seq._length; endPos++) {
			/*get the base*/
			ch = bases[endPos];

			/*unknown bases*/
			if (ch == BWT_NUM_NUCLEOTIDE) {
				lastStopAtN = true;
				range2 = make_uint2(1, 0);
				break;
			}

			//calculate the range
			range2.x = _rbwt->_bwtCCounts[ch] + _rbwt->bwtOcc(ch, range.x - 1)
					+ 1;
			range2.y = _rbwt->_bwtCCounts[ch] + _rbwt->bwtOcc(ch, range.y);
			if (range2.x > range2.y) {
				break;
			}
			range = range2;
			//Utils::log("%c", "ACGTN"[ch]);
		}
		/*Utils::log("\n startPos: %d endPos %d range: %u %u range2: %u %u\n", startPos, endPos, range.x, range.y, range2.x,
		 range2.y);*/

		/*If an exact match is found to the end of the query, no need to check any more*/
		if (range2.x <= range2.y) {
			//Utils::log("seed length %d min seed length %d #repeats: %d\n", endPos - startPos, minSeedSize, range2.y - range2.x + 1);
			if (range2.y - range2.x >= _maxSeedOcc) {
				range2.y = range2.x + _maxSeedOcc - 1;
			}
			ranges[numRanges++] = make_uint4(range2.x, range2.y, startPos,
					endPos - startPos);
			break;
		} else if (range.x <= range.y && endPos - startPos >= minSeedSize)/*a mismatch found and record the exact matches to the previous base*/
		{
			if (range.y - range.x >= _maxSeedOcc) {
				range.y = range.x + _maxSeedOcc - 1;
			}
			if (lastPos != endPos) {
				ranges[numRanges++] = make_uint4(range.x, range.y, startPos,
						endPos - startPos);
				/*update the last stop position*/
				lastPos = endPos;
				/*Utils::log("seed length %d startPos %d endPos %d min seed length %d #repeats: %d\n",
				 endPos - startPos, startPos, endPos, minSeedSize, range.y - range.x + 1);*/
			}

		}
		/*filter out consecutive Ns and re-locate the starting position*/
		if (lastStopAtN) {
			for (startPos = endPos + 1; startPos < numSuffixes; startPos++) {
				if (bases[startPos] != BWT_NUM_NUCLEOTIDE) {
					break;
				}
			}
		} else {
			startPos++;
		}
	}

	return numRanges;
}

uint32_t MemEngine::_genStripedSeeds(Sequence& seq, int strand, uint4* ranges,
		uint32_t seedSize, uint32_t stride) {
	uint8_t ch;
	uint2 range;
	uint32_t numRanges = 0;
	uint32_t startPos, endPos;
	uint8_t* bases = (strand == 0) ? seq._bases : seq._rbases;

	/*seeking seeds*/
	//Utils::log("strand: %d\n", strand);
	for (startPos = 0, endPos = seedSize; endPos < seq._length; startPos +=
			stride, endPos = startPos + seedSize) {

		/*for the current substring*/
		range = make_uint2(0, _bwtSeqLength);
		for (uint32_t pos = startPos; pos < endPos; ++pos) {
			/*get the base*/
			ch = bases[pos];
			//Utils::log("%c", "ACGTN"[ch]);

			/*unknown bases*/
			if (ch == BWT_NUM_NUCLEOTIDE) {
				range = make_uint2(1, 0);
				break;
			}

			//calculate the range
			range.x = _rbwt->_bwtCCounts[ch] + _rbwt->bwtOcc(ch, range.x - 1)
					+ 1;
			range.y = _rbwt->_bwtCCounts[ch] + _rbwt->bwtOcc(ch, range.y);
			if (range.x > range.y) {
				break;
			}
		}
		/*Utils::log("\nrange: %u %u %u\n", range.x, range.y,
		 range.y >= range.x ? range.y - range.x + 1 : 0);*/
		/*save the non-overlapping seeds*/
		if (range.x <= range.y) {
			if (range.y - range.x < 10 * _maxSeedOcc) {
				ranges[numRanges++] = make_uint4(range.x, range.y, startPos,
						seedSize);
			}
		}
	}

	return numRanges;
}

uint32_t MemEngine::_locateSAI(Sequence& seq, vector<uint4>& ranges,
		vector<uint4>& rranges, uint32_t& numRanges, uint32_t& numRranges,
		int seedType) {

	uint32_t minSeedSize;

	minSeedSize = _options->getMinSeedSize(seq._length);
	if (seq._tlength < minSeedSize) {
		return 0;
	}

	/*allocate space*/
	ranges.resize(seq._length);
	rranges.resize(seq._length);

	if (seedType == MAXIMAL_EXACT_MATCH_SEED) {
		/*for the forward strand*/
		numRanges = _genMEMSeeds(seq, 0, &ranges[0], minSeedSize);

		/*for the reverse strand*/
		numRranges = _genMEMSeeds(seq, 1, &rranges[0], minSeedSize);

		/*check the availability of seeds*/
		//Utils::log("line %d numRanges: %d numRranges: %d\n", __LINE__, numRanges, numRranges);
		if (numRanges + numRranges == 0) {
			/*for the forward strand*/
			numRanges = _genMEMSeeds(seq, 0, &ranges[0],
					(minSeedSize + GLOBAL_MIN_SEED_SIZE) >> 1);

			/*for the reverse strand*/
			numRranges = _genMEMSeeds(seq, 1, &rranges[0],
					(minSeedSize + GLOBAL_MIN_SEED_SIZE) >> 1);
		}
	} else {
		/*for the forward strand*/
		numRanges = _genStripedSeeds(seq, 0, &ranges[0], minSeedSize,
				minSeedSize);

		/*for the reverse strand*/
		numRranges = _genStripedSeeds(seq, 1, &rranges[0], minSeedSize,
				minSeedSize);
	}

	return numRanges + numRranges;
}
void MemEngine::_locateSeeds(Sequence& seq, vector<Seed>& seeds, int seedType) {
	Seed* seedsPtr;
	size_t numSeeds;
	int64_t target;
	uint4 *ranges, *rranges;
	uint32_t numFwSeeds, numRcSeeds;
	uint32_t numRanges, numRranges;

	/*clear all seeds*/
	seeds.clear();
	numFwSeeds = numRcSeeds = 0;
	/*calculate suffix array intervals*/
	if (_locateSAI(seq, _ranges, _rranges, numRanges, numRranges, seedType)
			== 0) {
		return;
	}

	/*calculate the total number of seeds*/
	/*for the forward strand*/
	ranges = &_ranges[0];
	for (uint4* p = ranges, *q = ranges + numRanges; p < q; ++p) {
		numFwSeeds += p->y - p->x + 1;
	}
	/*for the reverse strand*/
	rranges = &_rranges[0];
	for (uint4* p = rranges, *q = rranges + numRranges; p < q; ++p) {
		numRcSeeds += p->y - p->x + 1;
	}
	numSeeds = numFwSeeds + numRcSeeds;
	//Utils::log("numSeeds: %ld\n", numSeeds);

	/*allocate memory*/
	seeds.resize(numSeeds);
	seedsPtr = &seeds[0];

	/*calculate the mapping seeds for each occurrence of the seed*/
	Seed* t = seedsPtr;
	/*for the forward strand*/
	for (uint4 *p = ranges, *q = ranges + numRanges; p < q; ++p) {
		//for each suffix array interval
		for (uint32_t pos = p->x; pos <= p->y; ++pos) {
			/*get the mapping position in the reverse target sequence*/
			target = _rsa->getPosition(_rbwt, pos);
			/*get the mapping position in the forward target sequence*/
			target = _bwtSeqLength - target - p->w;

			/*save the mapping seeds*/
			t->_targetPosition = target;
			t->_queryPosition = p->z;
			t->_seedLength = p->w;
			t->_strand = 0;

			/*get genome index*/
			_rgenome->getGenomeIndex(target, t->_genomeIndex);

			/*increase the pointer address*/
			++t;
		}
	}
	/*for the reverse complement*/
	for (uint4 *p = rranges, *q = rranges + numRranges; p < q; ++p) {
		//for each suffix array interval
		for (uint32_t pos = p->x; pos <= p->y; ++pos) {
			/*get the mapping position in the reverse target sequence*/
			target = _rsa->getPosition(_rbwt, pos);
			/*get the mapping position in the forward target sequence*/
			target = _bwtSeqLength - target - p->w;

			/*save the mapping seeds*/
			t->_targetPosition = target;
			t->_queryPosition = p->z;
			t->_seedLength = p->w;
			t->_strand = 1;

			/*get genome index*/
			_rgenome->getGenomeIndex(target, t->_genomeIndex);

			/*increase the pointer address*/
			++t;
		}
	}

#if 0
	fprintf(stderr, "query, target, slen, strand, best genome\n");
	for (Seed* p = seedsPtr, *q = seedsPtr + numSeeds; p < q; ++p) {
		fprintf(stderr, "%u %u %u %u %d\n", p->_queryPosition, p->_targetPosition,
				p->_seedLength, p->_strand, p->_genomeIndex);
	}
	fprintf(stderr, "---------------#Seeds %ld -----------------------\n",
			numSeeds);

#endif
}

/*color space to base space conversion*/
pair<uint8_t*, int32_t> MemEngine::_cs2bsReadConvertBWA(const Sequence& seq,
		Mapping* mapping) {

	int length;
	uint8_t *newRead;
	uint8_t *csRead, *bsRead, *bsRef, *btArray, *tmpArray;
	CigarAlign * align;
	uint32_t* cigars;
	int32_t numCigars;

	if (!mapping) {
		return make_pair((uint8_t*) NULL, 0);
	}

	align = mapping->_align;
	cigars = align->getCigar();
	numCigars = align->getCigarLength();

	/*allocate intermediate buffers*/
	length = align->getAlignLength() + seq._length;
	_csBuffer.resize(length * 7);

	/*get the pointers*/
	csRead = &_csBuffer[0];
	bsRead = csRead + length;
	bsRef = bsRead + length;
	tmpArray = btArray = bsRef + length; /*share the same space*/

	/*extract sequences from CIGAR*/
	int32_t cigar, op;
	int64_t gPosition = mapping->_gposition; /*mapping position on the genome*/
	int basePosition = 0;
	int baseIndex = 0;
	int strand = mapping->_strand;

#if 0
	/*output the input cigar*/
	align->cigarOut(stderr);
	Utils::log("numMismathces: %d numGaps: %d numIndelds: %d\n",
			align->getNumMismatches(), align->getNumGaps(),
			align->getNumIndels());
	Utils::log("mapping->_gposition: %d\n", mapping->_gposition);
#endif

	/*get the colors*/
	uint8_t*bases = strand ? seq._rbases : seq._bases;
	/*set the leading base for the reference*/
	bsRef[0] = gPosition ? _getRefBase(_basePacGenome, gPosition) : 4;
	for (int32_t i = 0; i < numCigars; ++i) {
		cigar = cigars[i];
		length = cigar >> 2;
		op = cigar & 3;
		if (op == ALIGN_OP_M) {
			for (int j = 0; j < length; ++j) {
				/*get the color-space read base*/
				csRead[baseIndex] = _genCSBase(bases, seq._quals, strand,
						seq._length, basePosition);

				/*get the reference base*/
				bsRef[baseIndex + 1] = _getRefBase(_basePacGenome, gPosition);
				++baseIndex;

				/*move to the next base position*/
				++gPosition;
				++basePosition;
			}
		} else if (op == ALIGN_OP_I) {
			for (int j = 0; j < length; ++j) {
				/*get the color-space read base*/
				csRead[baseIndex] = _genCSBase(bases, seq._quals, strand,
						seq._length, basePosition);

				/*set to N or random base? FIXME*/
				bsRef[baseIndex + 1] = 4;
				++baseIndex;

				/*move to the next base position*/
				++basePosition;
			}
		} else if (op == ALIGN_OP_S) {
			/*soft clipping on the read*/
			basePosition += length;
		} else { /*op == ALIGN_OP_D*/
			gPosition += length;
		}
	}
	/*get the length of effective bases*/
	length = baseIndex;
//Utils::log("length: %d seqLength %d\n", length, seq._length);

	/*convert from color space to base space*/
	_convertCS2BS(csRead, bsRead, bsRef, btArray, length);
	newRead = _genBSQuals(bsRead, csRead, tmpArray, length);

	/*generate a new sequence*/
	uint8_t base, qual;
	uint8_t* data = new uint8_t[length * 3];
	if (!data) {
		Utils::exit("Memory allocation failed at line %d in file %s\n",
				__LINE__, __FILE__);
	}

	bases = strand ? data + length : data;
	uint8_t* rbases = strand ? data : data + length;
	uint8_t* quals = data + 2 * length;
	for (int32_t i = 0, j = length - 1; i < length; ++i, --j) {
		base = newRead[i] >> 6;
		qual = newRead[i] & 0x3f;
		if (qual == 63) {
			quals[strand ? j : i] = 33;
			bases[i] = rbases[j] = 4;
		} else {
			quals[strand ? j : i] = qual + 33;
			bases[i] = base;
			rbases[j] = _complements[base];
		}
	}

#if 0
	Utils::log("##alignment:");
	align->cigarOut(stderr);

	Utils::log("\n##Original color-space read:\n");
	for (int32_t i = 0; i < length; ++i) {
		if (i == 0) {
			Utils::log("%c", "ACGTN"[seq._bases[i]]);
		} else {
			Utils::log("%c", "0123."[seq._bases[i]]);
		}
	}
	Utils::log("\n##Translated base-space read:\n");
	int c1 = seq._bases[0];
	Utils::log("%c", "ACGTN"[c1]);
	for (int32_t i = 1; i < length; ++i) {
		int c2 = seq._bases[i];
		/*get the next base*/
		c2 = _cs2bsTable[c1][c2];
		Utils::log("%c", "ACGTN"[c2]);
		c1 = c2;
	}

	Utils::log("\n##BWA converted base-space read:\n");
	for (int32_t i = 0; i < length; ++i) {
		Utils::log("%c", "ACGTN"[data[i]]);
	}
	Utils::log("\n##Reference sequence:\n");
	for (int32_t i = 0; i < length; ++i) {
		Utils::log("%c", "ACGTN"[bsRef[i]]);
	}
	Utils::log("\n\n");
#endif

	return make_pair(data, length);
}

/*
 {A,C,G,T,N} -> {0,1,2,3,4}
 bsRef[0..size]: nucleotide reference: 0/1/2/3/4
 csRead[0..size-1]: color read+qual sequence: base<<6|qual; qual==63 for N
 bsRead[0..size]: nucleotide read sequence: 0/1/2/3 (returned)
 btArray[0..4*size]: backtrack array (working space)
 */
void MemEngine::_convertCS2BS(uint8_t* csRead, uint8_t* bsRead, uint8_t* bsRef,
		uint8_t* btArray, int32_t size) {
	int h[8], curr, last;
	int x, y, pos;
	int base, ref, qual, score;
	int minScore, minScoreBase;

// h[0..3] and h[4..7] are the current and last best score array, depending on curr and last
// recursion: initial value
	if (bsRef[0] > 3) {
		memset(h, 0, sizeof(int) << 2);
	} else {
		for (x = 0; x < 4; ++x) {
			h[x] = NUCL_MM;
		}
		h[bsRef[0]] = 0;
	}

// recursion: main loop
	curr = 1;
	last = 0;
	for (pos = 1; pos <= size; ++pos) {
		/*get the current base and its base quality score*/
		base = csRead[pos - 1] >> 6;
		qual = csRead[pos - 1] & 0x3f;

		/*get the current base on the reference*/
		ref = bsRef[pos];

		/*for each base at this position*/
		for (x = 0; x < 4; ++x) {

			/*select the best last base*/
			minScore = 0x7fffffff;
			minScoreBase = 0;
			for (y = 0; y < 4; ++y) {
				score = h[last << 2 | y];
				/*color mismatch penalty*/
				if (qual != 63 && base != _bs2csTable[1 << x | 1 << y]) {
					score += (qual < COLOR_MM) ? COLOR_MM : qual; // color mismatch
				}

				/*base mismatch penalty*/
				if (ref < 4 && ref != x) {
					score += NUCL_MM; // nt mismatch
				}

				/*select the best score*/
				if (score < minScore) {
					minScore = score;
					minScoreBase = y;
				}
			}
			/*save the best score for base x*/
			h[curr << 2 | x] = minScore;
			btArray[pos << 2 | x] = minScoreBase;
		}
		last = curr;
		curr = 1 - curr; // swap
	}

// back trace
	minScore = 0x7fffffff;
	minScoreBase = 0;
	for (base = 0; base < 4; ++base) {
		score = h[last << 2 | base];
		if (score < minScore) {
			minScore = score;
			minScoreBase = base;
		}
	}
	bsRead[size] = minScoreBase;

	for (pos = size - 1; pos >= 0; --pos) {
		bsRead[pos] = btArray[(pos + 1) << 2 | bsRead[pos + 1]];
	}
}

/*
 bsRead[0..size]: nucleotide read sequence: 0/1/2/3
 csRead[0..size-1]: color read+qual sequence: base<<6|qual; qual==63 for N
 tmpArray[0..size*2-1]: temporary array
 */
uint8_t* MemEngine::_genBSQuals(uint8_t *bsRead, uint8_t *csRead,
		uint8_t *tmpArray, int size) {
	int pos, c1, c2, qual, t1, t2;
	uint8_t *t2array = tmpArray + size;

	/*get the color sequence from converted bases*/
	c1 = bsRead[0];
	for (pos = 1; pos <= size; ++pos) {
		c2 = bsRead[pos]; // in principle, there is no 'N' in bsRead[]; just in case
		tmpArray[pos - 1] =
				(c1 >= 4 || c2 >= 4) ? 4 : _bs2csTable[1 << c1 | 1 << c2];
		c1 = c2;
	}

	/*starting from position 1*/
	for (pos = 1; pos <= size; ++pos) {
		qual = 0;
		c1 = csRead[pos - 1];
		c2 = csRead[pos];
		t1 = tmpArray[pos - 1];
		t2 = tmpArray[pos];
		if (t1 == (c1 >> 6) && t2 == (c2 >> 6)) {
			qual = (int) (c1 & 0x3f) + (int) (c2 & 0x3f) + 10;
		} else if (t1 == c1 >> 6) {
			qual = (int) (c1 & 0x3f) - (int) (c2 & 0x3f);
		} else if (t2 == c2 >> 6) {
			qual = (int) (c2 & 0x3f) - (int) (c1 & 0x3f);
		} // else, qual = 0
		qual = max(qual, 0);
		qual = min(qual, 60);

		/*save the converted base quality score*/
		t2array[pos] = bsRead[pos] << 6 | qual;
		if ((c1 & 0x3f) == 63 || (c2 & 0x3f) == 63) {
			t2array[pos] = 0;
		}
	}
	/*Keep the mapping quality score of the leading base*/
	t2array[0] = bsRead[0] << 6 | (csRead[0] & 0x3f);

	return t2array;
}

