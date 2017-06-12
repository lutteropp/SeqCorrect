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

#include <seqan/bam_io.h>
#include <seqan/sequence.h>
#include <string>
#include <vector>

#include "read_with_alignments.hpp"

namespace seq_correct {
namespace eval {

class BAMIterator {
public:
	BAMIterator();
	BAMIterator(const std::string &alignmentFilename);
	ReadWithAlignments next();
	std::vector<ReadWithAlignments> next(size_t numReads);
	bool hasReadsLeft();
	size_t numReadsLeft();
	double progress();
	size_t getNumReadsTotal();
	size_t getNumReadsTotalMapped();
	unsigned long long getNumReadsTotalUniqueMapped();
private:
	void countNumberOfReads(const std::string &alignmentFilename);
	seqan::BamFileIn bamFileIn;
	size_t readsLeft; size_t numReadsTotal; size_t numReadsTotalMapped; unsigned long long numReadsTotalUniqueMapped;

	std::vector<seqan::BamAlignmentRecord> records;seqan::CharString currentReadName;
};

} // end of namespace seq_correct::eval
} // end of namespace seq_correct
