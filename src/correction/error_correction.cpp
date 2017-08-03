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
#include "error_correction.hpp"
#include "../kmer/classification.hpp"

namespace seq_correct {
namespace correction {

using namespace classification;

size_t findSmallestNonrepetitive(const std::string& str, size_t pos, FMIndexMatcher& kmerCounter, PerfectUniformSequencingModel& pusm) {
	KmerType type = KmerType::REPEAT;
	for (size_t i = 1; i < str.size() - pos; i+=2) {
		type = classifyKmer(str.substr(pos, i), kmerCounter, pusm);
		if (type != KmerType::REPEAT) {
			return i;
		}
	}
	return std::numeric_limits<size_t>::max();
}

Read correctRead_kmer(const io::Read& read, FMIndexMatcher& kmerCounter, PerfectUniformSequencingModel& pusm,
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
Read correctRead_suffix_tree(const io::Read& read, FMIndexMatcher& kmerCounter, PerfectUniformSequencingModel& pusm,
		bool correctSingleIndels, bool correctMultidels) {
	throw std::runtime_error("not implemented yet");
}
Read correctRead_msa(const io::Read& read, FMIndexMatcher& kmerCounter, PerfectUniformSequencingModel& pusm,
		bool correctSingleIndels, bool correctMultidels) {
	throw std::runtime_error("not implemented yet");
}

} // end of namespace seq_correct::correction
} // end of namespace seq_correct
