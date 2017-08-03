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

#include <sdsl/csa_wt.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/rank_support_v5.hpp>
#include <sdsl/select_support_scan.hpp>
#include <sdsl/wt_huff.hpp>
#include <sdsl/suffix_arrays.hpp>

#include "../io/sequence_io.hpp"
#include "../external/const_string_ptr.hpp"

namespace seq_correct {
namespace counting {

using namespace sdsl;
typedef csa_wt<wt_huff<bit_vector, rank_support_v5<>, select_support_scan<>, select_support_scan<0>>, 1 << 20, 1 << 20> FMIndex;

class FMIndexMatcher {
public:
	FMIndexMatcher(const std::string& filename);
	size_t countKmer(const std::string &kmer);
	size_t countKmer(const external::ConstStringPtr& kmerPtr);
	size_t countKmerNoRC(const std::string &kmer);
	size_t countKmerNoRC(const external::ConstStringPtr& kmerPtr);
private:
	FMIndex fmIndex;
};

} // end of namespace seq_correct::counting
} // end of namespace seq_correct
