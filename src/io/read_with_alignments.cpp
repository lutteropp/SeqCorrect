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

#include "read_with_alignments.hpp"

#include <cassert>
#include <iostream>
#include <string>
#include <algorithm>
#include <stdexcept>

#include "../util/util.hpp"
#include "../util/enums.hpp"

namespace seq_correct {
namespace eval {

ReadWithAlignments::ReadWithAlignments(std::vector< seqan::BamAlignmentRecord> &bamRecords) :
		records(bamRecords) {
	// go through records, try to find one which has no hard-clipping
	bool found = false;
	for (seqan::BamAlignmentRecord record : bamRecords) {
		found = true;
		seqan::String<seqan::CigarElement<char> > cigar = record.cigar;
		for (size_t i = 0; i < length(cigar); ++i) {
			if (cigar[i].operation == 'H') {
				found = false;
				break;
			}
		}
		if (found) { // found a record which has no hard clipping
			name = toCString(record.qName);
			seq = "";
			qual = "";
			// copy the read
			for (size_t i = 0; i < length(record.seq); ++i) {
				seq += record.seq[i];
			}
			// if the read is reverse-complemented, reverse-complement it again to be normal.
			if (hasFlagRC(record)) {
				seq = util::reverseComplementString(seq);
				for (int i = length(record.qual) - 1; i >= 0; --i) {
					qual += record.qual[i];
				}
			} else {
				for (size_t i = 0; i < length(record.qual); ++i) {
					qual += record.qual[i];
				}
			}

			beginPos = record.beginPos;
			endPos = record.beginPos + seq.size() - 1;
			break;
		}
	}
	if (!found) {
		throw std::runtime_error("Read only occurs with hard clipping");
	}
}

} // end of namespace seq_correct::eval
} // end of namespace seq_correct
