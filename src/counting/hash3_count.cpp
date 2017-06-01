/**
 *  This file is part of pgrep.
 *
 *  pgrep is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  pgrep is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Copyright 2015-2016  Sarah Lutteropp, Bernhard Scheirle and Kevin Zerr
 */

#include <iostream>
#include <cassert>
#include "hash3_count.hpp"
#include "../util/util.hpp"

#define RANK3 3

namespace seq_correct {
namespace counting {

Hash3PreprocessInfo Hash3StringMatcher::preprocessPattern(const external::ConstStringPtr& pattern) {
	Hash3PreprocessInfo info;

	std::size_t m = (int) pattern.size();
	std::size_t mMinus1 = m - 1;
	std::size_t mMinus2 = m - 2;

	for (std::size_t a = 0; a < ALPHABET_SIZE; ++a) {
		info.shift[a] = mMinus2;
	}

	// note: do not try to write this hash code in a own function,
	// it will result in slower code even when marked with inline
	uint32_t h = pattern.getU(0);
	h = (h << 1) + pattern.getU(1);
	h = (h << 1) + pattern.getU(2);

	info.shift[h % ALPHABET_SIZE] = m - RANK3; // TODO: is the modulo neccessary?
	for (std::size_t i = RANK3; i < mMinus1; ++i) {
		h = pattern.getU(i - 2);
		h = (h << 1) + pattern.getU(i - 1);
		h = (h << 1) + pattern.getU(i);
		info.shift[h % ALPHABET_SIZE] = mMinus1 - i; // TODO: is the modulo neccessary?
	}

	h = pattern.getU(m - 3);
	h = (h << 1) + pattern.getU(mMinus2);
	h = (h << 1) + pattern.getU(mMinus1);

	info.shl = info.shift[h % ALPHABET_SIZE];
	if (info.shl == 0) {
		info.shl = 1;
	}
	info.shift[h % ALPHABET_SIZE] = 0;
	return info;
}

size_t Hash3StringMatcher::countString(const external::ConstStringPtr& pattern, const external::ConstStringPtr& string,
		const Hash3PreprocessInfo& info) {
	const size_t m = pattern.size();
	const size_t n = string.size();
	size_t count = 0;

	size_t j = m - 1;
	while (true) {
		size_t sh = 1;
		while (sh != 0 && j < n) {
			size_t h = string.getU(j - 2);
			h = (h << 1) + string.getU(j - 1);
			h = (h << 1) + string.getU(j);
			sh = info.shift[h % ALPHABET_SIZE];
			j += sh;
		}
		if (j < n) {
			assert(j + 1 >= m);

			size_t pos = j + 1 - m;
			size_t i = 0;
			while (i < m && static_cast<unsigned char>(pattern[i]) == string.getU(pos + i))
				i++;
			if (i >= m) {
				count++; // match found at position pos in the string
			}
			j += info.shl;
		} else {
			return count;
		}
	}
}

size_t Hash3StringMatcher::countInFile(const external::ConstStringPtr& pattern, const std::string& filepath,
		bool alsoReverseComplement) {
	size_t count = 0;
	Hash3PreprocessInfo infoPattern = preprocessPattern(pattern);
	io::ReadInput reader;
	reader.openFile(filepath);
	if (alsoReverseComplement) {
		std::string patternRC = util::reverseComplementString(pattern);
		external::ConstStringPtr patternRCPtr(&patternRC);
		Hash3PreprocessInfo infoPatternRC = preprocessPattern(patternRCPtr);
		while (reader.hasNext()) {
			io::Read read = reader.readNext(true, false, false);
			count += countString(pattern, external::ConstStringPtr(&read.seq), infoPattern);
			count += countString(patternRCPtr, external::ConstStringPtr(&read.seq), infoPatternRC);
		}
	} else {
		while (reader.hasNext()) {
			io::Read read = reader.readNext(true, false, false);
			count += countString(pattern, external::ConstStringPtr(&read.seq), infoPattern);
		}
	}
	return count;
}

size_t Hash3StringMatcher::countKmer(const std::string& kmer, const std::string& filepath) {
	external::ConstStringPtr kmerPtr(&kmer);
	return countInFile(kmerPtr, filepath, true);
}

size_t Hash3StringMatcher::countKmer(const external::ConstStringPtr& kmerPtr, const std::string& filepath) {
	return countInFile(kmerPtr, filepath, true);
}

size_t Hash3StringMatcher::countKmerNoRC(const std::string& kmer, const std::string& filepath) {
	external::ConstStringPtr kmerPtr(&kmer);
	return countInFile(kmerPtr, filepath, false);
}

size_t Hash3StringMatcher::countKmerNoRC(const external::ConstStringPtr& kmerPtr, const std::string& filepath) {
	return countInFile(kmerPtr, filepath, false);
}

} // end of namespace seq_correct::external
} // end of namespace seq_correct
