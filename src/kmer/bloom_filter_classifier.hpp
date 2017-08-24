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
#include "../external/bloom_filter.hpp"

namespace seq_correct {
namespace classification {

using namespace util;

/**
 * Class that keeps three bloom filters, one each for UNTRUSTED, UNIQUE, and REPETITIVE k-mers.
 * If a k-mer to-be-classified occurs in exactly *one* of the bloom filters, we instantly know its classification (with some little error).
 * If the k-mer occurs in multiple bloom filters, we need to reclassify the k-mer.
 * If the k-mer does not occur in any of the bloom filters, we classify the k-mer and add it to the corresponding filter.
 */

class BloomFilterClassifier {
public:
	BloomFilterClassifier(counting::Matcher& kmerCounter, pusm::PerfectUniformSequencingModel& pusm,
			coverage::CoverageBiasUnitMulti& biasUnit, const std::string& pathToOriginalReads);
	KmerType classifyKmer(const std::string& kmer);
private:
	counting::Matcher& kmerCounter;
	pusm::PerfectUniformSequencingModel& pusm;
	coverage::CoverageBiasUnitMulti& biasUnit;
	const std::string& pathToOriginalReads;

	bloom_filter filter_untrusted;
	bloom_filter filter_unique;
	bloom_filter filter_repeat;
};

inline BloomFilterClassifier::BloomFilterClassifier(counting::Matcher& kmerCounter,
		pusm::PerfectUniformSequencingModel& pusm, coverage::CoverageBiasUnitMulti& biasUnit,
		const std::string& pathToOriginalReads) :
		kmerCounter(kmerCounter), pusm(pusm), biasUnit(biasUnit), pathToOriginalReads(pathToOriginalReads) {
	// set up bloom filter for keeping track of already classified k-mers
	bloom_parameters parameters;
	parameters.projected_element_count = 1000000000; // expected number of k-mer to insert into the bloom filter
	parameters.false_positive_probability = 0.01;
	parameters.random_seed = 0xA5A5A5A5;
	parameters.compute_optimal_parameters();
	filter_untrusted = bloom_filter(parameters);
	filter_unique = bloom_filter(parameters);
	filter_repeat = bloom_filter(parameters);
}

inline KmerType BloomFilterClassifier::classifyKmer(const std::string& kmer) {
	bool inUntrusted = filter_untrusted.contains(kmer);
	bool inUnique = filter_unique.contains(kmer);
	bool inRepeat = filter_repeat.contains(kmer);
	uint8_t sum = inUntrusted + inUnique + inRepeat;
	if (sum == 1) {
		if (inUntrusted) {
			return KmerType::UNTRUSTED;
		} else if (inUnique) {
			return KmerType::UNIQUE;
		} else {
			return KmerType::REPEAT;
		}
	}

	KmerType type = classification::classifyKmer(kmer, kmerCounter, pusm, biasUnit, pathToOriginalReads);
	if (sum == 0) {
		if (type == KmerType::UNTRUSTED) {
			filter_untrusted.insert(kmer);
		} else if (type == KmerType::UNIQUE) {
			filter_unique.insert(kmer);
		} else {
			filter_repeat.insert(kmer);
		}
	}
	return type;
}


class BloomFilterClassifierGenome {
public:
	BloomFilterClassifierGenome(counting::Matcher& genomeIndex);
	KmerType classifyKmer(const std::string& kmer);
private:
	counting::Matcher& genomeIndex;

	bloom_filter filter_untrusted;
	bloom_filter filter_unique;
	bloom_filter filter_repeat;
};

inline BloomFilterClassifierGenome::BloomFilterClassifierGenome(counting::Matcher& genomeIndex) :
		genomeIndex(genomeIndex) {
	// set up bloom filter for keeping track of already classified k-mers
	bloom_parameters parameters;
	parameters.projected_element_count = 1000000000; // expected number of k-mer to insert into the bloom filter
	parameters.false_positive_probability = 0.01;
	parameters.random_seed = 0xA5A5A5A5;
	parameters.compute_optimal_parameters();
	filter_untrusted = bloom_filter(parameters);
	filter_unique = bloom_filter(parameters);
	filter_repeat = bloom_filter(parameters);
}

inline KmerType BloomFilterClassifierGenome::classifyKmer(const std::string& kmer) {
	bool inUntrusted = filter_untrusted.contains(kmer);
	bool inUnique = filter_unique.contains(kmer);
	bool inRepeat = filter_repeat.contains(kmer);
	uint8_t sum = inUntrusted + inUnique + inRepeat;
	if (sum == 1) {
		if (inUntrusted) {
			return KmerType::UNTRUSTED;
		} else if (inUnique) {
			return KmerType::UNIQUE;
		} else {
			return KmerType::REPEAT;
		}
	}

	uint16_t count = genomeIndex.countKmer(kmer);
	KmerType type = (count == 0) ? KmerType::UNTRUSTED : ((count == 1) ? KmerType::UNIQUE : KmerType::REPEAT);

	if (sum == 0) {
		if (type == KmerType::UNTRUSTED) {
			filter_untrusted.insert(kmer);
		} else if (type == KmerType::UNIQUE) {
			filter_unique.insert(kmer);
		} else {
			filter_repeat.insert(kmer);
		}
	}
	return type;
}

} // end of namespace seq_correct::kmer
} // end of namespace seq_correct
