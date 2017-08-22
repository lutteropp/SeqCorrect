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

#include "../pusm/pusm.hpp"
#include "../counting/fm_count.hpp"
#include "../util/enums.hpp"
#include "../coverage/coverage_bias.hpp"

namespace seq_correct {
namespace classification {

using namespace util;

/**
 * Classify a k-mer as either UNTRUSTED, UNIQUE, or REPETITIVE.
 * @param observedCount Count of the k-mer in the read dataset
 * @param expectedCount The expected count of the k-mer in an idealized sequencing setting
 * @param bias The G/C bias
 */
inline KmerType classifyKmer(size_t observedCount, double expectedCount, double bias) {
	double correctedCount = (double) 1.0 / bias * observedCount;
	if ((correctedCount < 2) || (correctedCount < 0.5 * expectedCount)) {
		return KmerType::UNTRUSTED;
	} else if (correctedCount <= 1.5 * expectedCount) {
		return KmerType::UNIQUE;
	} else {
		return KmerType::REPEAT;
	}
}

inline KmerType classifyKmerSarah(const std::string& kmer, size_t observedCount,
		pusm::PerfectUniformSequencingModel& pusm, coverage::CoverageBiasUnitSingle& biasUnit) {
	if (kmer.empty()) {
		throw std::runtime_error("The kmer is empty");
	}
	pusm::PusmData pusmData = pusm.expectedCount(kmer.size());
	return classifyKmer(observedCount, pusmData.expectation, biasUnit.computeCoverageBias(kmer));
}

inline KmerType classifyKmer(const std::string& kmer, counting::Matcher& kmerCounter,
		pusm::PerfectUniformSequencingModel& pusm, coverage::CoverageBiasUnitSingle& biasUnit) {
	return classifyKmerSarah(kmer, kmerCounter.countKmer(kmer), pusm, biasUnit);
}

inline KmerType classifyKmerReadBased(size_t k, size_t posInRead, const std::vector<uint16_t>& kmerCounts, double medianCount,
		const std::string& readSequence) {
	KmerType type;
	size_t count = kmerCounts[posInRead];
	if (count > medianCount * 1.5) {
		type = KmerType::REPEAT;
	} else if (count < medianCount * 0.5) {
		type = KmerType::UNTRUSTED;
	} else {
		type = KmerType::UNIQUE;
	}
	return type;
}

} // end of namespace seq_correct::kmer
} // end of namespace seq_correct
