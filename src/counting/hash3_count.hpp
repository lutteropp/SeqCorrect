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

#include <array>
#include "../external/const_string_ptr.hpp"
#include "matcher.hpp"

#pragma once

namespace seq_correct {
namespace counting {

#define ALPHABET_SIZE 5

struct Hash3PreprocessInfo {
	std::array<std::size_t, ALPHABET_SIZE> shift;
	std::size_t shl;
};

class Hash3StringMatcher : public Matcher {
public:
	Hash3StringMatcher(const std::string& filepath) : filepath(filepath) {};
	size_t countKmer(const std::string& kmer) override;
	size_t countKmer(const external::ConstStringPtr& kmerPtr) override;
	size_t countKmerNoRC(const std::string& kmer) override;
	size_t countKmerNoRC(const external::ConstStringPtr& kmerPtr) override;
private:
	std::string filepath;
	Hash3PreprocessInfo preprocessPattern(const external::ConstStringPtr& pattern);
	size_t countString(const external::ConstStringPtr& pattern, const external::ConstStringPtr& string, const Hash3PreprocessInfo& info);
	size_t countInFile(const external::ConstStringPtr& pattern, const std::string& filepath, bool alsoReverseComplement = true);
};

} // end of namespace seq_correct::counting
} // end of namespace seq_correct
