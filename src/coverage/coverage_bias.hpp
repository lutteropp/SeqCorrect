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
#include "../pusm/pusm.hpp"
#include "../io/sequence_io.hpp"
#include "../counting/fm_count.hpp"

namespace seq_correct {
namespace coverage {

/**
 * A struct holding gc-content and median coverage bias data.
 */
struct CoverageBiasData {
	double gc; //!< gc-content
	double bias; //!< median coverage bias for the given gc-content
};

/**
 * Computes the median coverage biases from the dataset for a specific k, then uses linear interpolation to respond to coverage bias queries.
 */
class CoverageBiasUnitSingle {
public:
	CoverageBiasUnitSingle() : gcStep(0) {};
	// preprocessing without a reference genome
	void preprocess(size_t k, const std::string &filepath, counting::Matcher& readsIndex, pusm::PerfectUniformSequencingModel &pusm);
	// preprocessing with a reference genome
	void preprocess(size_t k, const std::string& genome, counting::Matcher& readsIndex, counting::Matcher& genomeIndex);
	// using linear interpolation to compute coverage bias factor
	double computeCoverageBias(const std::string &kmer);
	double computeCoverageBias(double gc);
	void printMedianCoverageBiases();
private:
	std::vector<CoverageBiasData> medianCoverageBiases;
	double gcStep;
};

/**
 * Computes the median coverage biases from the dataset for multiple k, then uses linear interpolation to respond to coverage bias queries.
 */
class CoverageBiasUnitMulti {
public:
	CoverageBiasUnitMulti();
	// using linear interpolation to compute coverage bias factor
	double computeCoverageBias(const std::string &kmer, const std::string& filepath, counting::Matcher& readsIndex, pusm::PerfectUniformSequencingModel& pusm);
	double computeCoverageBias(const std::string &kmer, const std::string& genome, counting::Matcher& readsIndex, counting::Matcher& genomeIndex);
	void printMedianCoverageBiases();
private:
	std::unordered_map<size_t, CoverageBiasUnitSingle> biasUnits; //TODO: maybe change into an std::vector
};

} // end of namespace seq_correct::coverage
} // end of namespace seq_correct
