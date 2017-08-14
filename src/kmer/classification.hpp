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

namespace seq_correct {
namespace classification {

using namespace util;

/**
 * Classify a k-mer as either UNTRUSTED, UNIQUE, or REPETITIVE.
 * @param kmer The k-mer
 * @param pusmData The expected count and standard deviation of the k-mer in an idealized sequencing setting
 * @param observedCount Count of the k-mer in the read dataset
 */
inline KmerType classifyKmer(const pusm::PusmData& pusmData, size_t observedCount) {
	if (observedCount < 0.5 * pusmData.expectation) {
		return KmerType::UNTRUSTED;
	} else if (observedCount <= 1.5 * pusmData.expectation) {
		return KmerType::UNIQUE;
	} else {
		return KmerType::REPEAT;
	}
}

template<typename T>
inline KmerType classifyKmer(const T& kmer, counting::FMIndexMatcher& kmerCounter,
		pusm::PerfectUniformSequencingModel& pusm) {
	if (kmer.empty()) {
		throw std::runtime_error("The kmer is empty");
	}
	pusm::PusmData pusmData = pusm.expectedCount(kmer.size());
	size_t observedCount = kmerCounter.countKmer(kmer);
	return classifyKmer(pusmData, observedCount);
}

inline KmerType classifyKmerReadBased(size_t k, size_t posInRead, const std::vector<size_t>& kmerCounts, size_t medianCount,
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
