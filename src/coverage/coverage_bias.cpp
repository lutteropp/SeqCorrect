/*
 SeqCorrect - A toolkit for correcting Next Generation Sequencing data.
 Copyright (C) 2017 Sarah Lutteropp

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.

 Contact:
 Sarah Lutteropp <sarah.lutteropp@h-its.org>
 Exelixis Lab, Heidelberg Institute for Theoretical Studies
 Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
 */

#include <stdexcept>
#include <cmath>
#include <algorithm>
#include <unordered_set>
#include <iostream>
#include "coverage_bias.hpp"
#include "../external/bloom_filter.hpp"
#include "../util/util.hpp"
#include "../io/sequence_io.hpp"

namespace seq_correct {
namespace coverage {

#define USE_BLOOM_FILTER 1

/**
 * Fix biases with value 0 by applying linear interpolation. Zero-value biases occur if the dataset does not contain
 * any k-mers of a specific G/C-content.
 * @param data The median coverage bias values as learned from the read dataset
 */
void fixMissingBiases(std::vector<CoverageBiasData>& data) {
	auto nonZeroBias = [](const CoverageBiasData& dat) {return dat.bias != 0;};
	// find first index of nonzero bias, then fix entries data[0]...data[first-1]
	auto itFirst = std::find_if(data.begin(), data.end(), nonZeroBias);
	size_t first = std::distance(data.begin(), itFirst);
	if (first == data.size()) {
		for (size_t i = 0; i < data.size(); ++i) {
			data[i].bias = 1.0;
			std::cout << "All biases were missing. Fixed them by assuming no bias." << std::endl;
			return;
		}
	}
	if (first != 0) {
		for (size_t i = 0; i < first; ++i) {
			data[i].bias = data[first].bias;
		}
	}
	// find last index of nonzero bias, then fix entries data[last + 1]...data[data.size()-1]
	auto itLast = std::find_if(data.rbegin(), data.rend(), nonZeroBias);
	size_t last = std::distance(itLast, data.rend());
	if (last != data.size()) {
		for (size_t i = last; i < data.size(); ++i) {
			data[i].bias = data[last - 1].bias;
		}
	}
	// fix remaining possible zeros between data[first+1] and data[last-1]
	for (size_t i = first + 1; i < last; ++i) {
		if (data[i].bias == 0) {
			auto itLeft = std::find_if(data.rbegin(), data.rend() - i - 1, nonZeroBias);
			size_t left = std::distance(itLeft, data.rend()) - 1; // TODO: is the -1 correct here?
			auto itRight = std::find_if(data.begin() + i, data.end(), nonZeroBias);
			size_t right = std::distance(data.begin(), itRight);
			data[i].bias = data[left].bias + (data[right].bias - data[left].bias) / (data[right].gc - data[left].gc);
		}
	}
}

/**
 * Compute the median coverage biases, given the inferred biases from the dataset.
 * @param k the k-mer size
 * @param biases a vector containing vectors of bias factors as inferred from the dataset
 */
std::vector<CoverageBiasData> computeMedianBiases(size_t k, std::vector<std::vector<double> >& biases) {
	std::vector<CoverageBiasData> data(k+1);
	for (size_t i = 0; i < biases.size(); ++i) {
		std::sort(biases[i].begin(), biases[i].end());
		data[i].gc = i / (double) k;
		if (biases[i].empty()) {
			data[i].bias = 0.0;
		} else {
			size_t sizeHalved = biases[i].size() / 2;
			if (biases[i].size() % 2 == 0) {
				data[i].bias = (biases[i][sizeHalved] + biases[i][sizeHalved + 1]) / 2.0;
			} else {
				data[i].bias = biases[i][sizeHalved + 1];
			}
		}
	}
	fixMissingBiases(data);
	return data;
}

/**
 * Infer the coverage bias factor by comparing the count of a k-mer in the read dataset with the count of the k-mer
 * in the genome.
 * @param kmerPtr A pointer to the k-mer
 * @param readsIndex A structure for querying how often the k-mer occurs in the read dataset
 * @param genomeIndex A structure for querying how often the k-mer occurs in the genome
 */
double inferBias(const std::string& kmer, counting::Matcher& readsIndex,
		counting::Matcher& genomeIndex) {
	size_t countObserved = readsIndex.countKmer(kmer) + readsIndex.countKmer(util::reverseComplementString(kmer));
	size_t countGenome = genomeIndex.countKmer(kmer) + genomeIndex.countKmer(util::reverseComplementString(kmer));
	return countObserved / (double) countGenome;
}

/**
 * Infer the coverage bias factor by comparing the count of a k-mer in the read dataset with the expected count of
 * the k-mer in an idealized setting. Returns zero if the k-mer occurs very rarely in the dataset and is believed
 *  to be erroneous.
 * @param kmerPtr A pointer to the k-mer
 * @param readsIndex A structure for querying how often the k-mer occurs in the read dataset
 * @param pusm A structure for computing the expected count of k-mers in the read dataset, in an idealized
 * sequencing setting
 */
double inferBias(size_t k, const std::string& kmer, counting::Matcher& readsIndex,
		pusm::PerfectUniformSequencingModel &pusm) {
	size_t countObserved = readsIndex.countKmer(kmer);
	double countExpected = pusm.expectedCount(k).expectation;
	if (countObserved >= countExpected * 0.2) {
		// if this condition is left out, the coverage biases will be very low due to erroneous k-mers
		return countObserved / countExpected;
	} else {
		return 0.0;
	}
}

/**
 * Infer median coverage bias values from the read dataset.
 * Use a bloom filter to keep track of which k-mers have already been visited.
 * @param k The k-mer size
 * @param filepath The location of the file containing the reads
 * @param readsIndex A structure for querying how often the k-mer occurs in the read dataset
 * @param pusm A structure for computing the expected count of k-mers in the read dataset, in an idealized
 * sequencing setting
 */
std::vector<CoverageBiasData> preprocessWithoutGenome(size_t k, const std::string &filepath,
		counting::Matcher& readsIndex, pusm::PerfectUniformSequencingModel &pusm) {
	std::vector<std::vector<double> > biases;
	biases.resize(k + 1);

#ifdef USE_BLOOM_FILTER
	// set up bloom filter for keeping track of already visited k-mers
	bloom_parameters parameters;
	parameters.projected_element_count = 100000; // expected number of k-mer to insert into the bloom filter
	parameters.false_positive_probability = 0.01;
	parameters.random_seed = 0xA5A5A5A5;
	parameters.compute_optimal_parameters();
	bloom_filter filter(parameters);
#else
	std::unordered_set<std::string> visited;
#endif

	io::ReadInput reader;
	reader.openFile(filepath);

	while (reader.hasNext()) {
		io::Read seqRead = reader.readNext(true, false, false);
		std::string kmer = seqRead.seq.substr(0, k);
		size_t gcCount = util::countGC(kmer);
		double bias = inferBias(k, kmer, readsIndex, pusm);
		if (bias > 0) {
			biases[gcCount].push_back(bias);
		}
#ifdef USE_BLOOM_FILTER
		filter.insert(&seqRead.seq[0], k);
#else
		visited.insert(kmerPtr.toString());
#endif

		for (size_t i = 1; i + k < seqRead.seq.size(); ++i) {
			if (kmer[0] == 'G' || kmer[0] == 'C') {
				gcCount--;
			}
			kmer = seqRead.seq.substr(i, k);
			if (kmer[k - 1] == 'G' || kmer[k - 1] == 'C') {
				gcCount++;
			}

#ifdef USE_BLOOM_FILTER
			if (!filter.contains(&seqRead.seq[i], k)) {
#else
				if (visited.find(kmerPtr.toString()) == visited.end()) {
#endif
				double bias = inferBias(k, kmer, readsIndex, pusm);
				if (bias > 0) {
					biases[gcCount].push_back(bias);
				}
#ifdef USE_BLOOM_FILTER
				filter.insert(&seqRead.seq[i], k);
#else
				visited.insert(kmer);
#endif
			}
		}
	}
	return computeMedianBiases(k, biases);
}

/**
 * Infer median coverage bias values from the read dataset and the reference genome.
 * @param k The k-mer size
 * @param genome The reference genome
 * @param readsIndex A structure for querying how often the k-mer occurs in the read dataset
 * @param genomeIndex A structure for querying how often the k-mer occurs in the genome
 */
std::vector<CoverageBiasData> preprocessWithGenome(size_t k, const std::string& genome, counting::Matcher& readsIndex,
		counting::Matcher& genomeIndex) {
	std::vector<std::vector<double> > biases;
	biases.resize(k + 1);
	std::string kmer = genome.substr(0, k);
	size_t gcCount = util::countGC(kmer);
	double bias = inferBias(kmer, readsIndex, genomeIndex);
	if (bias > 0) {
		biases[gcCount].push_back(bias);
	}
	for (size_t i = 1; i + k < genome.size(); ++i) {
		if (kmer[0] == 'G' || kmer[0] == 'C') {
			gcCount--;
		}
		kmer = genome.substr(i, k);
		if (kmer[k - 1] == 'G' || kmer[k - 1] == 'C') {
			gcCount++;
		}
		double bias = inferBias(kmer, readsIndex, genomeIndex);
		if (bias > 0) {
			biases[gcCount].push_back(bias);
		}
	}
	return computeMedianBiases(k, biases);
}

/**
 * Infer median coverage biases from the dataset without using a reference genome.
 * @param k The kmer size used for training
 * @param filepath The location of the file containing the reads
 * @param readsIndex A structure for querying how often the k-mer occurs in the read dataset
 * @param pusm A structure for computing the expected count of k-mers in the read dataset, in an idealized
 * sequencing setting
 */
void CoverageBiasUnitSingle::preprocess(size_t k, const std::string &filepath, counting::Matcher& readsIndex,
		pusm::PerfectUniformSequencingModel& pusm) {
	gcStep = 1 / (double) k;
	medianCoverageBiases = preprocessWithoutGenome(k, filepath, readsIndex, pusm);
}

/**
 * Infer median coverage biases from the dataset by using a reference genome.
 * @param k The k-mer size used for training
 * @param genome The reference genome
 * @param readsIndex A structure for querying how often the k-mer occurs in the read dataset
 * @param genomeIndex A structure for querying how often the k-mer occurs in the genome
 */
void CoverageBiasUnitSingle::preprocess(size_t k, const std::string &genome, counting::Matcher& readsIndex,
		counting::Matcher& genomeIndex) {
	gcStep = 1 / (double) k;
	medianCoverageBiases = preprocessWithGenome(k, genome, readsIndex, genomeIndex);
}

/**
 * Computes the expected coverage bias of a k-mer, based on its GC-content.
 * @param gc The GC-content of the k-mer
 */
double CoverageBiasUnitSingle::computeCoverageBias(double gc) {
	size_t idxMin = std::floor(gc / gcStep);
	size_t idxMax = idxMin + 1;
	double gcMin = idxMin * gcStep;
	double gcMax = idxMax * gcStep;
	double biasMin = medianCoverageBiases[idxMin].bias;
	double biasMax = medianCoverageBiases[idxMax].bias;
	// linear interpolation, see https://en.wikipedia.org/wiki/Interpolation
	return biasMin + (biasMax - biasMin) / (gcMax - gcMin) * (gc - gcMin);
}

/**
 * Compute the expected coverage bias of a k-mer.
 * @param kmer The k-mer
 */
double CoverageBiasUnitSingle::computeCoverageBias(const std::string &kmer) {
	return computeCoverageBias(util::gcContent(kmer));
}

/**
 * Print the median coverage bias factors.
 */
void CoverageBiasUnitSingle::printMedianCoverageBiases() {
	for (size_t i = 0; i < medianCoverageBiases.size(); ++i) {
		std::cout << medianCoverageBiases[i].gc << " " << medianCoverageBiases[i].bias << "\n";
	}
}



CoverageBiasUnitMulti::CoverageBiasUnitMulti() {}

double CoverageBiasUnitMulti::computeCoverageBias(const std::string& kmer, const std::string& filepath, counting::Matcher& matcher, pusm::PerfectUniformSequencingModel& pusm) {
	size_t k = kmer.size();
	if (biasUnits.find(k) != biasUnits.end()) {
		return biasUnits[k].computeCoverageBias(kmer);
	} else {
		CoverageBiasUnitSingle cbsSingle;
		biasUnits[k] = cbsSingle;
		cbsSingle.preprocess(k, filepath, matcher, pusm);
		return cbsSingle.computeCoverageBias(kmer);
	}
}

double CoverageBiasUnitMulti::computeCoverageBias(const std::string &kmer, const std::string& genome, counting::Matcher& readsIndex, counting::Matcher& genomeIndex) {
	size_t k = kmer.size();
		if (biasUnits.find(k) != biasUnits.end()) {
			return biasUnits[k].computeCoverageBias(kmer);
		} else {
			CoverageBiasUnitSingle cbsSingle;
			biasUnits[k] = cbsSingle;
			cbsSingle.preprocess(k, genome, readsIndex, genomeIndex);
			return cbsSingle.computeCoverageBias(kmer);
		}
}

void CoverageBiasUnitMulti::printMedianCoverageBiases() {
	for (auto& kv : biasUnits) {
		std::cout << kv.first << ":\n";
		kv.second.printMedianCoverageBiases();
		std::cout << "\n";
	}
}


} // end of namespace seq_correct::coverage
} // end of namespace seq_correct
