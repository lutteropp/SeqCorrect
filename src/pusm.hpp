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
#include "info.hpp"

/**
 * Perfect Uniform Sequencing Model (PUSM). Computes the expected count and standard deviation of a DNA-sequence of
 * length k in the read dataset. For computing these values, we need the estimated genome size and the read length
 * distribution.
 *
 * The PUSM assumes the following conditions:
 * - the k-mer is unique in the genome
 * - there is no coverage bias
 * - the reads map perfectly to the genome, without any errors or non-genomic content
 */

namespace pusm {

// TODO: Or is it better to just use an std::pair<double, double>?
// TODO: Is this the right way to go if I want to avoid confusing a standard deviation with an expected count?
class PusmData {
public:
	PusmData(double stdev, double expectation) :
			_stdev(stdev), _expectation(expectation) {
	}
	double stdev();
	double expectation();
private:
	double _stdev;
	double _expectation;
};

// TODO: where to put this one? Global variable seems bad...
std::unordered_map<size_t, PusmData> pusmBuffer; // buffer storing already computed answers for given k-mer lengths

PusmData expectedCount(info::GenomeType type, size_t genomeSize, const std::unordered_map<size_t, size_t> &readLengths, size_t k);

} // end of namespace pusm
