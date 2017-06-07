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

#pragma once

#include <unordered_map>
#include <algorithm>
#include "../util/util.hpp"
#include "../util/genome_type.hpp"

namespace seq_correct {
namespace pusm {

/**
 * A struct that stores standard deviation and expected count of a k-mer in the read dataset.
 */
struct PusmData {
	PusmData() :
			stdev(0), expectation(0) {
	}
	double stdev;
	double expectation;
};

/**
 * Perfect Uniform Sequencing Model (PUSM). Computes the expected count and standard deviation of a DNA-sequence of
 * length k in the read dataset. For computing these values, we need the estimated genome size and the read length
 * distribution.
 *
 * The PUSM assumes the following conditions:
 * - the k-mer is unique in the genome
 * - there is no coverage bias
 * - the reads map perfectly to the genome, without any errors or non-genomic content
 * - the reads contain no adapter sequences or other artificial sequences
 */
class PerfectUniformSequencingModel {
public:
	PerfectUniformSequencingModel(util::GenomeType type, size_t genomeSize,
			const std::unordered_map<size_t, size_t>& readLengths) :
			type(type), genomeSize(genomeSize), readLengths(readLengths) {
	}
	PusmData expectedCount(size_t k);
private:
	util::GenomeType type;
	size_t genomeSize;
	std::unordered_map<size_t, size_t> readLengths;
	std::unordered_map<size_t, PusmData> pusmBuffer; // buffer storing already computed answers for given k-mer lengths
};

// TODO: Decide on whether these functions should be hidden (this is, only defined in the .cpp file) or not
PusmData expectedCountLinear(size_t genomeSize, const std::unordered_map<size_t, size_t>& readLengths, size_t k);
PusmData expectedCountCircular(size_t genomeSize, const std::unordered_map<size_t, size_t>& readLengths, size_t k);

} // end of namespace seq_correct::pusm
} // end of namespace seq_correct
