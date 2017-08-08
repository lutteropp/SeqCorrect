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

#include <stdexcept>
#include <functional>
#include "error_correction.hpp"
#include "../kmer/classification.hpp"

namespace seq_correct {
namespace correction {

using namespace classification;
using namespace std::placeholders;

Read correctRead_none(const Read& read, FMIndexMatcher& kmerCounter, PerfectUniformSequencingModel& pusm,
		bool correctSingleIndels = true, bool correctMultidels = false);
Read correctRead_simple_kmer(const Read& read, FMIndexMatcher& kmerCounter, PerfectUniformSequencingModel& pusm,
		bool correctSingleIndels = true, bool correctMultidels = false);
Read correctRead_adaptive_kmer(const Read& read, FMIndexMatcher& kmerCounter, PerfectUniformSequencingModel& pusm,
		bool correctSingleIndels = true, bool correctMultidels = false);
Read correctRead_suffix_tree(const io::Read& read, FMIndexMatcher& kmerCounter, PerfectUniformSequencingModel& pusm,
		bool correctSingleIndels = true, bool correctMultidels = false);
Read correctRead_full_msa(const Read& read, FMIndexMatcher& kmerCounter, PerfectUniformSequencingModel& pusm,
		bool correctSingleIndels = true, bool correctMultidels = false);
Read correctRead_partial_msa(const Read& read, FMIndexMatcher& kmerCounter, PerfectUniformSequencingModel& pusm,
		bool correctSingleIndels = true, bool correctMultidels = false);

size_t findSmallestNonrepetitive(const std::string& str, size_t pos, FMIndexMatcher& kmerCounter,
		PerfectUniformSequencingModel& pusm) {
	KmerType type = KmerType::REPEAT;
	for (size_t i = 1; i < str.size() - pos; i += 2) {
		type = classifyKmer(str.substr(pos, i), kmerCounter, pusm);
		if (type != KmerType::REPEAT) {
			return i;
		}
	}
	return std::numeric_limits<size_t>::max();
}

Read correctRead_none(const io::Read& read, FMIndexMatcher& kmerCounter, PerfectUniformSequencingModel& pusm,
		bool correctSingleIndels, bool correctMultidels) {
	(void) kmerCounter;
	(void) pusm;
	(void) correctSingleIndels;
	(void) correctMultidels;
	return read;
}

Read correctRead_simple_kmer(const io::Read& read, FMIndexMatcher& kmerCounter, PerfectUniformSequencingModel& pusm,
		bool correctSingleIndels, bool correctMultidels) {
	/*
	 * TODO:
	 * - use sliding window of fixed size? Or just... increase k-mer size as long as k-mer is repetitive?
	 * - check whether highly similar k-mers are also trusted/ other errors are also highly supported
	 * --> if so, further increase k-mer size
	 * - find a good way to deal with conflicting correction candidates from different iterations (this was missing in thesis)
	 * - find a better minimum size for k than 1
	 */

	io::Read correctedRead(read);
	size_t pos = 0;
	while (pos < correctedRead.seq.size()) { // loop over starting position
		//find smallest nonrepetitive k-mer size
		size_t k = findSmallestNonrepetitive(correctedRead.seq, pos, kmerCounter, pusm);
		if (k == std::numeric_limits<size_t>::max()) {
			break;
		}
		KmerType type = classifyKmer(correctedRead.seq.substr(pos, k), kmerCounter, pusm);
		if (type == KmerType::UNIQUE)

			pos++;
	}

	throw std::runtime_error("not implemented yet");
}

Read correctRead_adaptive_kmer(const io::Read& read, FMIndexMatcher& kmerCounter, PerfectUniformSequencingModel& pusm,
		bool correctSingleIndels, bool correctMultidels) {
	throw std::runtime_error("not implemented yet");
}

Read correctRead_suffix_tree(const io::Read& read, FMIndexMatcher& kmerCounter, PerfectUniformSequencingModel& pusm,
		bool correctSingleIndels, bool correctMultidels) {
	throw std::runtime_error("not implemented yet");
}

Read correctRead_full_msa(const io::Read& read, FMIndexMatcher& kmerCounter, PerfectUniformSequencingModel& pusm,
		bool correctSingleIndels, bool correctMultidels) {
	throw std::runtime_error("not implemented yet");
}

Read correctRead_partial_msa(const io::Read& read, FMIndexMatcher& kmerCounter, PerfectUniformSequencingModel& pusm,
		bool correctSingleIndels, bool correctMultidels) {
	throw std::runtime_error("not implemented yet");
}

void correctReads(const std::string& readsFilepath, CorrectionAlgorithm algo, FMIndexMatcher& kmerCounter,
		PerfectUniformSequencingModel& pusm, const std::string& outputPath) {

	bool correctSingleIndels = true;
	bool correctMultidels = false;

	std::function<Read(const Read&, FMIndexMatcher&, PerfectUniformSequencingModel&, bool, bool)> f;
	switch (algo) {
	case CorrectionAlgorithm::NONE:
		f = std::bind(correctRead_none, _1, kmerCounter, pusm, correctSingleIndels, correctMultidels);
		break;
	case CorrectionAlgorithm::SIMPLE_KMER:
		f = std::bind(correctRead_simple_kmer, _1, kmerCounter, pusm, correctSingleIndels, correctMultidels);
		break;
	case CorrectionAlgorithm::ADAPTIVE_KMER:
		f = std::bind(correctRead_adaptive_kmer, _1, kmerCounter, pusm, correctSingleIndels, correctMultidels);
		break;
	case CorrectionAlgorithm::FULL_MSA:
		f = std::bind(correctRead_full_msa, _1, kmerCounter, pusm, correctSingleIndels, correctMultidels);
		break;
	case CorrectionAlgorithm::PARTIAL_MSA:
		f = std::bind(correctRead_partial_msa, _1, kmerCounter, pusm, correctSingleIndels, correctMultidels);
		break;
	case CorrectionAlgorithm::SUFFIX_TREE:
		f = std::bind(correctRead_suffix_tree, _1, kmerCounter, pusm, correctSingleIndels, correctMultidels);
		break;
	default:
		throw std::runtime_error("Unknown correction algorithm");
	}

	io::ReadOutput printer;
	printer.createFile(outputPath);

	io::ReadInput reader;
	reader.openFile(readsFilepath);

#pragma omp parallel
	{
#pragma omp single
		{
			while (reader.hasNext()) {
				io::Read uncorrected = reader.readNext(true, true, true);
#pragma omp task shared(printer, kmerCounter, pusm, correctSingleIndels, correctMultidels) firstprivate(uncorrected)
				{
					io::Read corrected = f(uncorrected, kmerCounter, pusm, correctSingleIndels, correctMultidels);
#pragma omp critical
					printer.write(corrected);
				}
			}
		}
#pragma omp taskwait
	}
}

} // end of namespace seq_correct::correction
} // end of namespace seq_correct
