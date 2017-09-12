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

#include <unordered_map>

#include "../io/sequence_io.hpp"
#include "enums.hpp"

namespace seq_correct {
namespace util {

inline char reverseComplementChar(char c) {
	if (c == 'A') {
		return 'T';
	} else if (c == 'C') {
		return 'G';
	} else if (c == 'G') {
		return 'C';
	} else if (c == 'T') {
		return 'A';
	} else {
		return c;
	}
}

template<typename T>
inline std::string reverseComplementString(const T& text) {
	std::string revComp = "";
	for (size_t i = 0; i < text.size(); ++i) {
		revComp += reverseComplementChar(text[text.size() - i - 1]);
	}
	return revComp;
}

io::Read reverseComplementRead(const io::Read& read);
double gcContent(const std::string& kmer);

template<typename T>
inline size_t countGC(const T& kmer) {
	size_t gcCount = 0;
	for (size_t i = 0; i < kmer.size(); ++i) {
		if (kmer[i] == 'G' || kmer[i] == 'C') {
			gcCount++;
		}
	}
	return gcCount;
}

inline bool isSelfReverseComplement(const std::string& kmer) {
	if (kmer.size() % 2 == 1) {
		return false;
	}
	bool same = true;
	for (size_t i = 0; i < kmer.size() / 2; ++i) {
		if (kmer[i] != util::reverseComplementChar(kmer[kmer.size() - i - 1])) {
			same = false;
			break;
		}
	}
	return same;
}

inline std::unordered_map<size_t, size_t> countReadLengths(const std::string& readFilepath) {
	std::unordered_map<size_t, size_t> res;
	io::ReadInput reader(readFilepath);
	while (reader.hasNext()) {
		io::Read read = reader.readNext(true, false, false);
		if (res.find(read.seq.size()) == res.end()) {
			res[read.seq.size()] = 1;
		} else {
			res[read.seq.size()]++;
		}
	}
	return res;
}

std::string kmerAfterError(const std::string& kmer, size_t pos, ErrorType type);

size_t kmerToNumber(const std::string& kmer);
std::string numberToKmer(size_t n, size_t k);

struct pairhash {
public:
  template <typename T, typename U>
  std::size_t operator()(const std::pair<T, U> &x) const
  {
    return 13 * std::hash<T>()(x.first) + 37 * std::hash<U>()(x.second);
  }
};

template <typename T, typename U>
inline std::ostream& operator<<(std::ostream& os, const std::pair<T, U> &x)
{
	return os << x.first << ' ' << x.second;
}

template <typename T, typename U>
inline std::istream& operator>>(std::istream& is, std::pair<T, U> &x)
{
	return is >> x.first >> x.second;
}


} // end of namespace seq_correct::util
} // end of namespace seq_correct
