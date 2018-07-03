/*
 * PairedEnd.cpp
 *
 *  Created on: Jan 18, 2012
 *      Author: yongchao
 */

#include "PairedEnd.h"
#include "SeqFileParser.h"
#include "Options.h"
#include <omp.h>

PairedEnd::PairedEnd(Options* options, Genome * genome, SAM* sam) {
	_options = options;
	_genome = genome;
	_sam = sam;

	/*get the number of threads*/
	_numThreads = _options->getNumThreads();

	/*create parameters for threads*/
	_threads.resize(_numThreads);
	for (size_t tid = 0; tid < _threads.size(); ++tid) {
		_threads[tid] = new ThreadParams(tid, _options, _sam,
				new MemEngine(_options, _genome, _sam), this);
	}
	pthread_mutex_init(&_mutex, NULL);
}

PairedEnd::~PairedEnd() {
	pthread_mutex_destroy(&_mutex);
	for (size_t i = 0; i < _threads.size(); ++i) {
		delete _threads[i];
	}
	_threads.clear();
}
void PairedEnd::estimateInsertSize(int minAlignedPairs, int maxReadBatchSize,
		int numReadBatchs) {
	size_t nreads;
	int mapQualReliable = _options->getMapQualReliable();
	int pairingMode = _options->getPairingMode();
	int numReadPairs = 0;
	SeqFileParser *parser1, *parser2;
	Sequence *seqs1, *seqs2;
	vector<int> globalInsertSizes;
	int* localInsertSizes;
	double stime = Utils::getSysTime();

	Utils::log("Start estimating insert size using the top 0x%x read pairs\n",
			maxReadBatchSize * numReadBatchs);
	/*allocate memory*/
	seqs1 = new Sequence[maxReadBatchSize];
	if (!seqs1) {
		Utils::exit("Memory allocation failed in line %d in function %s\n",
				__LINE__, __FUNCTION__);
	}
	seqs2 = new Sequence[maxReadBatchSize];
	if (!seqs2) {
		Utils::exit("Memory allocation failed in line %d in function %s\n",
				__LINE__, __FUNCTION__);
	}
	localInsertSizes = new int[maxReadBatchSize];
	if (!localInsertSizes) {
		Utils::exit("Memory allocation failed in line %d in function %s\n",
				__LINE__, __FUNCTION__);
	}
	/*reserved space for global insert sizes*/
	globalInsertSizes.reserve(maxReadBatchSize * numReadBatchs);

	/*set the number of threads for OpenMP runtime*/
	omp_set_num_threads(_numThreads);

	/*for paired-end alignment*/
	bool done = false;
	vector<pair<string, int> > &inputs = _options->getInputFileList();
	for (size_t file = 0; file < inputs.size(); file += 2) {
		/*open the file for the left sequences*/
		parser1 = new SeqFileParser(_options, inputs[file].first.c_str(), false,
				inputs[file].second);
		/*open the file for the right sequences*/
		parser2 = new SeqFileParser(_options, inputs[file + 1].first.c_str(),
				false, inputs[file].second);

		/*read a batch of paired-end reads*/
		if ((nreads = parser1->getSeqLockFree(seqs1, maxReadBatchSize)) == 0) {
			Utils::log("Empty input file\n");
		}
		if (parser2->getSeqLockFree(seqs2, maxReadBatchSize) != nreads) {
			Utils::exit("The two files have different number of reads\n");
		}

		/*start the main loop*/
		do {
			size_t index;
			/*get the single-end alignments*/
#pragma omp parallel for private(index) default(shared) schedule(dynamic, 1)
			for (index = 0; index < nreads; ++index) {

				/*get the thread ID*/
				int tid = omp_get_thread_num();

				/*get the engine for the thread*/
				MemEngine *engine = _threads[tid]->_engine;

				/*perform alignment*/
				Mapping *mapping1, *mapping2;
				engine->align(seqs1[index], mapping1);
				engine->align(seqs2[index], mapping2);

				/*calcualte the insert size*/
				int insertSize = -1; /*dummy insert size*/
				if (mapping1 && mapping2) {
					if ((pairingMode == MATE_PAIRED_READS ?
							mapping1->_strand == mapping2->_strand :
							mapping1->_strand != mapping2->_strand)
							&& mapping1->_genomeIndex == mapping2->_genomeIndex
							&& mapping1->_mapQual >= mapQualReliable
							&& mapping2->_mapQual >= mapQualReliable) {
						if (pairingMode == MATE_PAIRED_READS) {
							insertSize = ((int64_t) mapping1->_position)
									- mapping2->_position;
						} else {
							if (mapping1->_strand == 0) {
								insertSize = ((int64_t) mapping1->_position)
										- mapping2->_position
										- seqs2[index]._length;
							} else {
								insertSize = ((int64_t) mapping1->_position)
										- mapping2->_position
										+ seqs1[index]._length;
							}
						}
						if (insertSize < 0) {
							insertSize = -insertSize;
						}
					}
				}
				/*save the insert size for the current read pair*/
				localInsertSizes[index] = insertSize;
				/*release the mapping*/
				if (mapping1)
					delete mapping1;
				if (mapping2)
					delete mapping2;
			}

			/*merge all local insert sizes*/
			for (size_t i = 0; i < nreads; ++i) {
				if (localInsertSizes[i] > 0) {
					globalInsertSizes.push_back(localInsertSizes[i]);
				}
			}

			/*statistical information*/
			numReadPairs += nreads;
			if (numReadPairs >= maxReadBatchSize * numReadBatchs) {
				Utils::log("#read pairs read from the input: %d\n",
						numReadPairs);
				done = true;
				break;
			}

			/*re-load a batch of read pairs*/
			if ((nreads = parser1->getSeqLockFree(seqs1, maxReadBatchSize))
					== 0) {
				break;
			}
			if (parser2->getSeqLockFree(seqs2, maxReadBatchSize) != nreads) {
				Utils::exit("The two files have different number of reads\n");
			}
		} while (1);

		/*release the file parser*/
		delete parser1;
		delete parser2;

		/*check if sufficient aligned pairs have got*/
		if (done) {
			break;
		}
	}
	/*release resources*/
	delete[] localInsertSizes;
	delete[] seqs2;
	delete[] seqs1;

	/*check the number of insert sizes*/
	int numInsertSizes = globalInsertSizes.size();
	if (numInsertSizes < minAlignedPairs) {
		Utils::exit(
				"#qualified reads pairs (%d) is less than %d and please specify the insert size through parameters\n",
				numInsertSizes, minAlignedPairs);
	}

	/*sort the insert sizes*/
	sort(globalInsertSizes.begin(), globalInsertSizes.end());

	/*get the inset size for percentile 25, 50 and 75*/
	int* insertSizes = &globalInsertSizes[0];
	int p25 = insertSizes[(int) (numInsertSizes * 0.25 + 0.499)];
	int p50 = insertSizes[(int) (numInsertSizes * 0.5 + 0.499)];
	int p75 = insertSizes[(int) (numInsertSizes * 0.75 + 0.499)];

	/*estimate the mean using mean value instead of the median? FIXME*/
	double dvariance = 0;
	int64_t mean, variance;
	int low, high, count;
	low = p25 - 2 * (p75 - p25);
	if (low < 0)
		low = 0;

	high = p75 + 2 * (p75 - p25);
	if (high > 5 * p50)
		high = 5 * p50;

	/*calculate the mean value*/
	mean = 0;
	count = 0;
	for (int i = 0; i < numInsertSizes; ++i) {
		if (insertSizes[i] > low && insertSizes[i] < high) {
			mean += insertSizes[i];
			++count;
		}
	}
	if (count > 0) {
		mean = (int64_t) (static_cast<double>(mean) / count + 0.499);
	} else {
		Utils::exit(
				"Failed to estimate the insert size. Please specify this information through parameters\n");
	}

	/*calculate the variance*/
	variance = 0;
	count = 0;
	for (int i = 0; i < numInsertSizes; ++i) {
		if (insertSizes[i] > low && insertSizes[i] < high) {
			variance += (insertSizes[i] - mean) * (insertSizes[i] - mean);
			++count;
		}
	}
	if (count > 0) {
		dvariance = sqrt(static_cast<double>(variance) / count);
	}
	if (variance == 0)
		dvariance = 1;

	/*further limit the variance*/
	if (dvariance > 0.2 * mean) {
		dvariance = 0.2 * mean;
	}
	variance = static_cast<int>(dvariance + 0.499);

	/*set the insert size as well as the standard deviation*/
	_options->setInsertSize(mean);
	_options->setStdInsertSize(variance);
	Utils::log("Estimated insert size: %ld +/- %ld from %d effective samples\n",
			mean, variance, count);

	/*release resource*/
	globalInsertSizes.clear();

	/*update distance information*/
	for (size_t tid = 0; tid < _threads.size(); ++tid) {
		_threads[tid]->_engine->updateDistance();
	}

	double etime = Utils::getSysTime();
	Utils::log(
			"Finish estimating insert size (taken %f seconds using %d threads)\n",
			etime - stime, _numThreads);
}
void PairedEnd::execute() {
	SeqFileParser *parser1, *parser2;
	vector<pthread_t> threadIDs;

	/*estimate the insert size?*/
	if (_options->estimateInsertSize()) {
		estimateInsertSize(100, INS_SIZE_EST_MULTIPLE,
				_options->getTopReadsEstIns() / INS_SIZE_EST_MULTIPLE);
	}

	/*initialize variables*/
	threadIDs.resize(_numThreads);

	/*for paired-end alignment*/
	vector<pair<string, int> > &inputs = _options->getInputFileList();
	for (size_t file = 0; file < inputs.size(); file += 2) {
		/*open the file for the left sequences*/
		parser1 = new SeqFileParser(_options, inputs[file].first.c_str(), false,
				inputs[file].second);
		/*open the file for the right sequences*/
		parser2 = new SeqFileParser(_options, inputs[file + 1].first.c_str(),
				false, inputs[file].second);

		/*create threads*/
		for (size_t tid = 0; tid < _threads.size(); ++tid) {
			/*set file parser*/
			_threads[tid]->setFileParser(parser1, parser2);

			/*create the thread*/
			if (pthread_create(&threadIDs[tid], NULL, _threadFunc,
					_threads[tid]) != 0) {
				Utils::exit("Thread creating failed\n");
			}
		}
		/*wait for the completion of all threads*/
		for (size_t tid = 0; tid < _threads.size(); ++tid) {
			pthread_join(threadIDs[tid], NULL);
		}
		/*release the file parser*/
		delete parser1;
		delete parser2;
	}
	threadIDs.clear();

	/*report the alignment information*/
	int64_t sum = 0, sum2 = 0, sum3 = 0;
	for (size_t tid = 0; tid < _threads.size(); ++tid) {
		sum += _threads[tid]->_numReads;
		sum2 += _threads[tid]->_numAligned;
		sum3 += _threads[tid]->_numPaired;
	}
	Utils::log("#Reads aligned: %ld / %ld (%.2f%%)\n", sum2, sum,
			((double) (sum2)) / sum * 100);

	Utils::log("#Reads paired: %ld / %ld (%.2f%%)\n", sum3, sum,
			((double) (sum3)) / sum * 100);
}
void* PairedEnd::_threadFunc(void* arg) {
	Sequence seq1, seq2;
	Mapping *mapping1, *mapping2;
	int64_t numReads = 0, numAligned = 0, numPaired = 0;
	ThreadParams *params = (ThreadParams*) arg;
	SAM* sam = params->_sam;
	Options* options = params->_options;
	SeqFileParser *parser1 = params->_parser1;
	SeqFileParser *parser2 = params->_parser2;
	PairedEnd* aligner = (PairedEnd*) params->_aligner;
	MemEngine *engine = params->_engine;
	double stime = Utils::getSysTime();
	double etime;
	bool aligned;
	int multi = options->getMaxMultiAligns();
	int minMapQual = options->getMinMapQual();
	vector<Mapping*> mapv1, mapv2;

	/*reserve space for the vectors*/
	mapv1.reserve(multi);
	mapv2.reserve(multi);
	bool paired, mapped1, mapped2;
	int flags1 = SAM_FPD | SAM_FR1; /*paired-end reads and the first read*/
	int flags2 = SAM_FPD | SAM_FR2; /*paired-end reads and the second read*/

	while (1) {
		/*read a sequence pair*/
		aligner->lock();
		if (!parser1->getSeqLockFree(seq1)) {
			aligner->unlock();
			break;
		}
		if (!parser2->getSeqLockFree(seq2)) {
			Utils::log("The two files have different number of sequences\n");
			aligner->unlock();
			break;
		}
		aligner->unlock();

		/*invoke the engine to get the paired-end alignments*/
		if ((paired = engine->align(seq1, seq2, mapv1, mapv2))) {
			numPaired++;
		}
		/*output the alignment*/
		options->lock();
		if (paired) { /*paired alignments found*/
			if (mapv1.size() != mapv2.size()) {
				Utils::exit(
						"Error occured for paired-end alignments in line %d in funciton %s\n",
						__LINE__, __FUNCTION__);
			}
			aligned = false;
			for (size_t i = 0; i < mapv1.size(); ++i) {
				mapping1 = mapv1[i];
				mapping2 = mapv2[i];
				if (mapping1->_mapQual >= minMapQual
						|| mapping2->_mapQual >= minMapQual) {
					sam->printPaired(seq1, *mapping1, *mapping2, true, flags1);
					sam->printPaired(seq2, *mapping2, *mapping1, true, flags2);
					aligned = true;
				} else {
					/*deemed to be unaligned*/
					sam->print(seq1, flags1 | SAM_FMU);
					sam->print(seq2, flags2 | SAM_FMU);
				}
			}
			if (aligned) {
				numAligned += 2;
			}
		} else {
			if (mapv1.size() >= 2 || mapv2.size() >= 2) {
				Utils::exit(
						"Error occured for paired-end alignments in line %d in funciton %s\n",
						__LINE__, __FUNCTION__);
			}
			/*get the mappings*/
			mapping1 = mapv1.size() > 0 ? mapv1[0] : NULL;
			mapping2 = mapv2.size() > 0 ? mapv2[0] : NULL;

			/*print out the mappings*/
			mapped1 = mapping1 && mapping1->_mapQual >= minMapQual;
			mapped2 = mapping2 && mapping2->_mapQual >= minMapQual;
			if (mapped1 && mapped2) {
				sam->printPaired(seq1, *mapping1, *mapping2, false, flags1);
				sam->printPaired(seq2, *mapping2, *mapping1, false, flags2);
				numAligned += 2;
			} else {
				if (mapped1) {
					sam->print(seq1, *mapping1, flags1 | SAM_FMU);
					++numAligned;
				} else {
					/*deemed to be unaligned*/
					if (mapped2) {
						sam->printPaired(seq1, *mapping2, flags1);
					} else {
						sam->print(seq1, flags1 | SAM_FMU);
					}
				}
				if (mapped2) {
					sam->print(seq2, *mapping2, flags2 | SAM_FMU);
					++numAligned;
				} else {
					if (mapped1) {
						sam->printPaired(seq2, *mapping1, flags2);
					} else {
						sam->print(seq2, flags2 | SAM_FMU);
					}
				}
			}
		}
		options->unlock();

		/*release the resources*/
		for (size_t i = 0; i < mapv1.size(); ++i) {
			delete mapv1[i];
		}
		mapv1.clear();
		for (size_t i = 0; i < mapv2.size(); ++i) {
			delete mapv2[i];
		}
		mapv2.clear();

		/*statistical information*/
		numReads++;
		if (numReads % 100000 == 0) {
			etime = Utils::getSysTime();
			Utils::log(
					"processed %ld read pairs by thread %d in %.2f seconds\n",
					numReads, params->_tid, etime - stime);
		}
	}
	/*return the alignment results*/
	params->_numAligned += numAligned;
	params->_numReads += numReads * 2;
	params->_numPaired += numPaired * 2;

	return NULL;
}

