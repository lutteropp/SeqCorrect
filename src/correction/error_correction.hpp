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

#include <vector>

#include "../util/util.hpp"
#include "../io/sequence_io.hpp"
#include "../counting/matcher.hpp"
#include "../pusm/pusm.hpp"

namespace seq_correct {
namespace correction {

using namespace io;
using namespace counting;
using namespace pusm;

Read correctRead_kmer(const Read& read, Matcher& kmerCounter, PerfectUniformSequencingModel& pusm, bool correctSingleIndels = true,
		bool correctMultidels = false);
Read correctRead_suffix_tree(const io::Read& read, Matcher& kmerCounter, PerfectUniformSequencingModel& pusm, bool correctSingleIndels = true, bool correctMultidels = false);
Read correctRead_msa(const Read& read, Matcher& kmerCounter, PerfectUniformSequencingModel& pusm, bool correctSingleIndels = true, bool correctMultidels = false);

} // end of namespace seq_correct::correction
} // end of namespace seq_correct
