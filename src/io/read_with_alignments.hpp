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
#include <string>
#include <seqan/bam_io.h>
#include "../util/enums.hpp"

/*
 * TODO: Bei einem Chimeric break dürfte die beginPos des Reads unterschiedlich sein!!!
 */

namespace seq_correct {
namespace eval {

using namespace util;

class ReadWithAlignments {
public:
	ReadWithAlignments(std::vector<seqan::BamAlignmentRecord> &bamRecords);
	std::string name;
	std::string seq;
	std::string qual;
	std::vector<seqan::BamAlignmentRecord> records;
	size_t beginPos;
	size_t endPos;
};

} // end of namespace seq_correct::eval
} // end of namespace seq_correct
