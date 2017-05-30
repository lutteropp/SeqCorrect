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

#include <vector>
#include "pusm.hpp"
#include "sequence_io.hpp"

namespace seq_correct {
namespace coverage {

class CoverageBiasData {
public:
	CoverageBiasData(double gc, double bias) :
			_gc(gc), _bias(bias) {
	}
	double gc() {
		return _gc;
	}
	double bias() {
		return _bias;
	}
private:
	double _gc;
	double _bias;
};

class CoverageBiasUnit {
public:
	CoverageBiasUnit();
	void preprocess(const std::vector<io::Read> &reads, pusm::PerfectUniformSequencingModel &pusm);
	void preprocess(const std::vector<io::Read> &reads, const std::string &genome);
	void preprocess(const std::vector<std::string> &filepaths, pusm::PerfectUniformSequencingModel &pusm);
	void preprocess(const std::vector<std::string> &filepaths, const std::string &genome);
	void preprocess(const std::string &filepath, pusm::PerfectUniformSequencingModel &pusm);
	void preprocess(const std::string &filepath, const std::string &genome);

	double computeCoverageBias(const std::string &kmer);
	double computeCoverageBias(size_t k, double gc);
private:
	std::vector<CoverageBiasData> _medianCoverageBiases;
};

} // end of namespace seq_correct::coverage
} // end of namespace seq_correct
