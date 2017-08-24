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

#include "classification.hpp"
#include "../external/sparsepp/spp.h"

namespace seq_correct {
namespace classification {

using namespace util;
using spp::sparse_hash_map;

/**
 * Class that buffers already classified k-mers.
 */

class HashClassifier {
public:
	HashClassifier(counting::Matcher& kmerCounter, pusm::PerfectUniformSequencingModel& pusm,
			coverage::CoverageBiasUnitMulti& biasUnit, const std::string& pathToOriginalReads);
	KmerType classifyKmer(const std::string& kmer);
private:
	counting::Matcher& kmerCounter;
	pusm::PerfectUniformSequencingModel& pusm;
	coverage::CoverageBiasUnitMulti& biasUnit;
	const std::string& pathToOriginalReads;
	sparse_hash_map<std::string, KmerType> buffer;
};

inline HashClassifier::HashClassifier(counting::Matcher& kmerCounter, pusm::PerfectUniformSequencingModel& pusm,
		coverage::CoverageBiasUnitMulti& biasUnit, const std::string& pathToOriginalReads) :
		kmerCounter(kmerCounter), pusm(pusm), biasUnit(biasUnit), pathToOriginalReads(pathToOriginalReads) {
}

inline KmerType HashClassifier::classifyKmer(const std::string& kmer) {
	if (buffer.contains(kmer)) {
		return buffer[kmer];
	}
	KmerType type = classification::classifyKmer(kmer, kmerCounter, pusm, biasUnit, pathToOriginalReads);
#pragma omp critical
	buffer[kmer] = type;
	return type;
}

class HashClassifierGenome {
public:
	HashClassifierGenome(counting::Matcher& genomeIndex);
	KmerType classifyKmer(const std::string& kmer);
private:
	counting::Matcher& genomeIndex;
	sparse_hash_map<std::string, KmerType> buffer;
};

inline HashClassifierGenome::HashClassifierGenome(counting::Matcher& genomeIndex) :
		genomeIndex(genomeIndex) {}

inline KmerType HashClassifierGenome::classifyKmer(const std::string& kmer) {
	if (buffer.contains(kmer)) {
		return buffer[kmer];
	}
	uint16_t count = genomeIndex.countKmer(kmer);
	KmerType type = (count == 0) ? KmerType::UNTRUSTED : ((count == 1) ? KmerType::UNIQUE : KmerType::REPEAT);
#pragma omp critical
	buffer[kmer] = type;
	return type;
}

} // end of namespace seq_correct::kmer
} // end of namespace seq_correct
