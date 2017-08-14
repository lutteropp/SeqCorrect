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
#include <string>
#include "../util/enums.hpp"
#include "../util/correction.hpp"
#include "../io/read_with_alignments.hpp"
#include "../counting/fm_count.hpp"
#include "../pusm/pusm.hpp"
#include "../coverage/coverage_bias.hpp"
#include "error_evaluation_data.hpp"

namespace seq_correct {
namespace eval {

struct HandlingInfo {
	unsigned positionInRead = 0;
	unsigned hardClippedBases = 0;
	unsigned softClippedBases = 0;
	unsigned insertedBases = 0;
	unsigned deletedBases = 0;

	bool revComp;
	unsigned beginPos;
	unsigned readLength;
	std::vector<Correction> corrections;
};

std::vector<Correction> extractErrors(ReadWithAlignments& rwa, const std::string &referenceGenome);
std::vector<Correction> extractErrors(const std::string& correctedRead, const std::string& referenceGenome);
std::vector<Correction> extractErrors(const std::string& correctedRead, const std::string& referenceGenome,
		size_t beginPos);

void updateEvaluationData(ErrorEvaluationData& data, const std::vector<Correction>& errorsTruth,
		const std::vector<Correction>& errorsPredicted, size_t readLength);

// Compare the corrected reads with the aligned original reads
/*EvaluationData evaluateCorrections(const std::string& originalReadsFilepath, const std::string& correctedReadsFilepath,
 const std::string& genomeFilepath);*/

ErrorEvaluationData evaluateCorrectionsByAlignment(const std::string& alignmentFilepath,
		const std::string& correctedReadsFilepath, const std::string& genomeFilepath);

KmerEvaluationData classifyKmersTestSarah(size_t k, util::GenomeType genomeType, const std::string& alignmentFilepath,
		const std::string& pathToOriginalReads, const std::string& genomeFilepath, counting::Matcher& fmReads, counting::Matcher& fmGenome);

KmerEvaluationData classifyKmersTestReadbased(size_t k, GenomeType genomeType, const std::string& alignmentFilepath, const std::string& pathToOriginalReads,
		const std::string& genomeFilepath, counting::Matcher& fmReads, counting::Matcher& fmGenome);

} // end of namespace seq_correct::eval
} // end of namespace seq_correct
