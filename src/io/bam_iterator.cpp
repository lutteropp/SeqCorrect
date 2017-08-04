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

#include <seqan/bam_io.h>
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <iostream>
#include <stdexcept>

#include "bam_iterator.hpp"
#include "read_with_alignments.hpp"

namespace seq_correct {
namespace eval {

void BAMIterator::countNumberOfReads(const std::string &alignmentFilename) {
	numReadsTotal = 0;
	numReadsTotalMapped = 0;
	seqan::BamFileIn bamFileInReadCounting;
	if (!seqan::open(bamFileInReadCounting, alignmentFilename.c_str())) {
		throw std::runtime_error("ERROR: Could not open " + alignmentFilename);
	}
	// read and discard header
	seqan::BamHeader header;
	seqan::readHeader(header, bamFileInReadCounting);
	seqan::BamAlignmentRecord record;
	if (!atEnd(bamFileInReadCounting)) {
		seqan::readRecord(record, bamFileInReadCounting);
		currentReadName = record.qName;
		if (!hasFlagUnmapped(record)) {
			numReadsTotalMapped++;
		}
		numReadsTotal++;
	}
	int countUnique = 1;
	while (!atEnd(bamFileInReadCounting)) {
		seqan::readRecord(record, bamFileInReadCounting);
		if (record.qName != currentReadName) {
			if (countUnique == 1) {
				numReadsTotalUniqueMapped++;
			}
			countUnique = 1;
			numReadsTotal++;
			if (!hasFlagUnmapped(record)) {
				numReadsTotalMapped++;
			}
			currentReadName = record.qName;
		} else {
			countUnique++;
		}
	}

	if (countUnique == 1) {
		numReadsTotalUniqueMapped++;
	}

	std::cout << "Finished counting reads. There are " << numReadsTotal << " reads in total.\n";
}

BAMIterator::BAMIterator() {
	numReadsTotal = 0;
	numReadsTotalMapped = 0;
	numReadsTotalUniqueMapped = 0;
	readsLeft = 0;
}

BAMIterator::BAMIterator(const std::string &alignmentFilename) {
	countNumberOfReads(alignmentFilename);
	readsLeft = numReadsTotal;

	// Open input file, BamFileIn can read SAM and BAM files.
	if (!seqan::open(bamFileIn, alignmentFilename.c_str())) {
		throw std::runtime_error("ERROR: Could not open " + alignmentFilename);
	}

	// read and discard header
	seqan::BamHeader header;
	seqan::readHeader(header, bamFileIn);
	seqan::BamAlignmentRecord record;

	if (!atEnd(bamFileIn)) {
		seqan::readRecord(record, bamFileIn);
		currentReadName = record.qName;
		records.push_back(record);
	}
}

bool BAMIterator::hasReadsLeft() {
	return (readsLeft > 0);
}

size_t BAMIterator::numReadsLeft() {
	return readsLeft;
}

double BAMIterator::progress() {
	return ((numReadsTotal - readsLeft) / ((double) numReadsTotal)) * 100;
}

// Assumes that the BAM file is sorted by read name
ReadWithAlignments BAMIterator::next() {
	//std::cout << "next called.\n";
	seqan::BamAlignmentRecord record;

	while (!atEnd(bamFileIn)) {
		seqan::readRecord(record, bamFileIn);
		if (record.qName == currentReadName) {
			records.push_back(record);
		} else {
			ReadWithAlignments alignedRead(records);
			readsLeft--;
			currentReadName = record.qName;

			records.clear();
			records.push_back(record);
			return alignedRead;
		}
	}

	if (records.empty()) {
		throw std::runtime_error("The records are empty!");
	}

	ReadWithAlignments alignedRead(records);
	readsLeft--;
	records.clear();

	assert(readsLeft == 0);
	return alignedRead;
}

std::vector<ReadWithAlignments> BAMIterator::next(size_t numReads) {
	std::vector<ReadWithAlignments> res;
	while (res.size() < numReads && readsLeft > 0) {
		res.push_back(next());
	}
	return res;
}

size_t BAMIterator::getNumReadsTotal() {
	return numReadsTotal;
}

size_t BAMIterator::getNumReadsTotalMapped() {
	return numReadsTotalMapped;
}

unsigned long long BAMIterator::getNumReadsTotalUniqueMapped() {
	return numReadsTotalUniqueMapped;
}

} // end of namespace seq_correct::eval
} // end of namespace seq_correct
