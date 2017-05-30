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
#include <iostream>
#include "coverage_bias.hpp"
#include "../external/constStringPtr.hpp"
#include "../counting/count_kmer.hpp"
#include "../util/util.hpp"

namespace seq_correct {
namespace coverage {

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
	if (last != data.size() - 1) {
		for (size_t i = last + 1; i < data.size(); ++i) {
			data[i].bias = data[last].bias;
		}
	}
	// fix remaining possible zeros between data[first+1] and data[last-1]
	for (size_t i = first + 1; i < last; ++i) {
		if (data[i].bias == 0) {
			auto itLeft = std::find_if(data.rbegin(), data.rend() - i - 1, nonZeroBias);
			size_t left = std::distance(itLeft, data.rend());
			auto itRight = std::find_if(data.begin() + i, data.end(), nonZeroBias);
			size_t right = std::distance(data.begin(), itRight);
			data[i].bias = data[left].bias + (data[right].bias - data[left].bias) / (data[right].gc - data[left].gc);
		}
	}
}

double inferBias(const external::ConstStringPtr& kmerPtr, counting::FMIndex& readsIndex,
		counting::FMIndex& genomeIndex) {
	size_t countObserved = readsIndex.countKmer(kmerPtr);
	size_t countGenome = genomeIndex.countKmer(kmerPtr);
	return countObserved / (double) countGenome;
}

std::vector<CoverageBiasData> preprocessWithoutGenome(size_t k, pusm::PerfectUniformSequencingModel &pusm) {
	throw std::runtime_error("not implemented yet");
}

std::vector<CoverageBiasData> preprocessWithGenome(size_t k, const std::string& genome, counting::FMIndex& readsIndex,
		counting::FMIndex& genomeIndex) {
	std::vector<CoverageBiasData> data;
	std::vector<std::vector<double> > biases;
	biases.resize(k + 1);
	external::ConstStringPtr genomePtr(&genome);
	external::ConstStringPtr kmerPtr = genomePtr.substr(0, k);
	size_t gcCount = util::countGC(kmerPtr);
	for (size_t i = 1; i + k < genome.size(); ++i) {
		if (kmerPtr[0] == 'G' || kmerPtr[0] == 'C') {
			gcCount--;
		}
		kmerPtr = genomePtr.substr(i, k);
		if (kmerPtr[k - 1] == 'G' || kmerPtr[k - 1] == 'C') {
			gcCount++;
		}
		biases[gcCount].push_back(inferBias(kmerPtr, readsIndex, genomeIndex));
	}
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

void CoverageBiasUnit::preprocess(size_t k, const std::string &filepath, counting::FMIndex& readsIndex,
		pusm::PerfectUniformSequencingModel& pusm) {
	_medianCoverageBiases = preprocessWithoutGenome(k, pusm);
}

void CoverageBiasUnit::preprocess(size_t k, const std::string &genome, counting::FMIndex& readsIndex,
		counting::FMIndex& genomeIndex) {
	_medianCoverageBiases = preprocessWithGenome(k, genome, readsIndex, genomeIndex);
}

double CoverageBiasUnit::computeCoverageBias(double gc) {
	size_t idxMin = std::floor(gc / _gcStep);
	size_t idxMax = idxMin + 1;
	double gcMin = idxMin * _gcStep;
	double gcMax = idxMax * _gcStep;
	double biasMin = _medianCoverageBiases[idxMin].bias;
	double biasMax = _medianCoverageBiases[idxMax].bias;
	// linear interpolation, see https://en.wikipedia.org/wiki/Interpolation
	return biasMin + (biasMax - biasMin) / (gcMax - gcMin) * (gc - gcMin);
}

double CoverageBiasUnit::computeCoverageBias(const std::string &kmer) {
	return computeCoverageBias(util::gcContent(kmer));
}

} // end of namespace seq_correct::coverage
} // end of namespace seq_correct
