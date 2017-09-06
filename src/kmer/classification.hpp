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
	expectedCount *= 2; // because of the reverse complement handlings

	double correctedCount = (double) 1.0 / bias * observedCount;
	if ((correctedCount < 2) || (correctedCount < 0.5 * expectedCount)) {
		return KmerType::UNTRUSTED;
	} else if (correctedCount <= 1.5 * expectedCount) {
		return KmerType::UNIQUE;
	} else {
		return KmerType::REPEAT;
	}
}

inline KmerType classifyKmer(const std::string& kmer, counting::Matcher& readsIndex,
		pusm::PerfectUniformSequencingModel& pusm, coverage::CoverageBiasUnitMulti& biasUnit,
		const std::string& pathToOriginalReads) {

	if (classifyKmer(readsIndex.countKmer(kmer), pusm.expectedCount(kmer.size()).expectation,
			biasUnit.computeCoverageBias(kmer, pathToOriginalReads, readsIndex, pusm)) == KmerType::REPEAT && kmer.size() > 50) {
		std::cout << readsIndex.countKmer(kmer) << "\n" << pusm.expectedCount(kmer.size()).expectation << "\n" << biasUnit.computeCoverageBias(kmer, pathToOriginalReads, readsIndex, pusm) << "\n";
		throw std::runtime_error("BGF");
	}

	return classifyKmer(readsIndex.countKmer(kmer), pusm.expectedCount(kmer.size()).expectation,
			biasUnit.computeCoverageBias(kmer, pathToOriginalReads, readsIndex, pusm));
}

inline KmerType classifyKmerReadBased(size_t k, size_t posInRead, const std::vector<uint16_t>& kmerCounts,
		double medianCount, const std::string& readSequence) {
	return classifyKmer(kmerCounts[posInRead], medianCount, 1.0);
}

} // end of namespace seq_correct::kmer
} // end of namespace seq_correct
