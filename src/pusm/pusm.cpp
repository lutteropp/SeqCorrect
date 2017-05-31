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

#include "pusm.hpp"
#include <stdexcept>
#include <cmath>

namespace seq_correct {
namespace pusm {

/**
 * Compute the expected count of a unique k-mer under the Perfect Uniform Sequencing Model with a linear genome.
 * @param genomeSize The (estimated) size of the genome
 * @param readLengths The read length distribution
 * @param k The k-mer size
 */
PusmData expectedCountLinear(size_t genomeSize, const std::unordered_map<size_t, size_t>& readLengths, size_t k) {
	double expected = 0.0;
	double variance = 0.0;
	for (auto pair : readLengths) {
		double l = pair.first;
		if (l > genomeSize) {
			throw std::runtime_error("l = " + std::to_string(l) + " > " + std::to_string(genomeSize));
		}
		double n = pair.second;
		if (l < k)
			continue;
		// had to be fixed, because in my Master thesis I forgot that the interval [a,b] contains
		//  b-a+1 elements and not b-a elements. -.-
		double p = ((l - k) * (genomeSize - l) + 1) / ((genomeSize - k + 1) * (genomeSize - l + 1));
		if (p > 1.0 || p != p) {
			throw std::runtime_error("Something is wrong with the probability");
		};
		expected += n * p;
		variance += n * p * (1 - p);
	}
	PusmData res;
	res.expectation = expected;
	res.stdev = sqrt(variance);
	return res;
}

/**
 * Compute the expected count of a unique k-mer under the Perfect Uniform Sequencing Model with a circular genome.
 * @param genomeSize The (estimated) size of the genome
 * @param readLengths The read length distribution
 * @param k The k-mer size
 */
PusmData expectedCountCircular(size_t genomeSize, const std::unordered_map<size_t, size_t>& readLengths, size_t k) {
	double expected = 0.0;
	double variance = 0.0;
	for (auto pair : readLengths) {
		double l = pair.first;
		if (l > genomeSize) {
			throw std::runtime_error("l = " + std::to_string(l) + " > " + std::to_string(genomeSize));
		}
		double n = pair.second;
		if (l < k)
			continue;
		// This formula is correct. The probability is NOT 1.0 in the case of read length = genome length,
		//  because the reads are not circular/ on a torus.
		double p = (double) (l - k + 1) / genomeSize;
		if (p > 1.0 || p != p) {
			throw std::runtime_error("Something is wrong with the probability");
		};
		expected += n * p;
		variance += n * p * (1 - p);
	}
	PusmData res;
	res.expectation = expected;
	res.stdev = sqrt(variance);
	return res;
}

/**
 * Compute the expected count of a unique k-mer under the Perfect Uniform Sequencing model.
 * @param k The k-mer size
 */
PusmData PerfectUniformSequencingModel::expectedCount(size_t k) {
	if (_pusmBuffer.find(k) != _pusmBuffer.end()) {
		return _pusmBuffer[k];
	}
	PusmData res;
	if (_type == util::GenomeType::CIRCULAR) {
		res = expectedCountCircular(_genomeSize, _readLengths, k);
	} else {
		res = expectedCountLinear(_genomeSize, _readLengths, k);
	}
	_pusmBuffer[k] = res;
	return res;
}

} // end of namespace seq_correct::pusm
} // end of namespace seq_correct
