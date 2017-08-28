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
#include "../counting/fm_count.hpp"
#include "../pusm/pusm.hpp"
#include "../util/enums.hpp"
#include "../coverage/coverage_bias.hpp"

namespace seq_correct {
namespace correction {

using namespace io;
using namespace counting;
using namespace pusm;
using namespace util;

std::string findReplacement(const std::string& kmer, ErrorType errorType, size_t posInKmer);

std::string kmerAfterError(const std::string& kmer, ErrorType error, int posOfError);

uint8_t numUntrustedKmers(const std::string& read, size_t minK, Matcher& kmerCounter, PerfectUniformSequencingModel& pusm,
		coverage::CoverageBiasUnitMulti& biasUnit, const std::string& pathToOriginalReads);

std::vector<uint8_t> badKmerCoverage(const std::string& read, size_t minK, Matcher& kmerCounter, PerfectUniformSequencingModel& pusm,
		coverage::CoverageBiasUnitMulti& biasUnit, const std::string& pathToOriginalReads);

bool readIsPerfect(const std::string& read, size_t minK, Matcher& kmerCounter, PerfectUniformSequencingModel& pusm,
		coverage::CoverageBiasUnitMulti& biasUnit, const std::string& pathToOriginalReads);

void correctReads(const std::string& pathToOriginalReads, CorrectionAlgorithm algo, Matcher& kmerCounter,
		PerfectUniformSequencingModel& pusm, coverage::CoverageBiasUnitMulti& biasUnit, const std::string& outputPath);

std::pair<size_t, size_t> affectedReadArea(size_t posInRead, size_t readLength, size_t kmerSize);

} // end of namespace seq_correct::correction
} // end of namespace seq_correct
