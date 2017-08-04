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

#include <stdexcept>
#include <fstream>
#include <algorithm>
#include <seqan/bam_io.h>
#include "evaluation.hpp"
#include "../util/util.hpp"
#include "../io/bam_iterator.hpp"
#include "../io/read_with_alignments.hpp"
#include "../io/sequence_io.hpp"
#include "alignment.hpp"

namespace seq_correct {
namespace eval {

std::string extractCigarString(const seqan::String< seqan::CigarElement<char> >& cigar) {
	std::string cigarString = "";
	for (unsigned i = 0; i < length(cigar); ++i) {
		cigarString += cigar[i].operation;
		cigarString += ":";
		cigarString += std::to_string(cigar[i].count);
		cigarString += " ";
	}
	return cigarString;
}

void handleChimericBreak(ReadWithAlignments& rwa, size_t cigarCount, HandlingInfo& info) {
	/*unsigned nucleotidePositionInReference = info.beginPos + info.positionInRead - info.insertedBases
	 - info.softClippedBases + info.deletedBases;*/
	unsigned realPositionInRead = info.positionInRead + info.hardClippedBases; // the position of the first base of the chimeric break

	assert(realPositionInRead < readLength);

	if (realPositionInRead == 0) { // in this case, take the position of the last base of the chimeric break
		realPositionInRead += cigarCount;
	}

	// What if multiple chimeric breaks happen in a read? ... It should be fine, too.
	//  As long as there aren't multiple hard clipped regions in a single record.
	info.hardClippedBases += cigarCount;
	if (info.revComp) {
		realPositionInRead = info.readLength - realPositionInRead - 1;
	}

	// because we treat chimeric breaks as deletion of multiple bases
	info.corrections.push_back(Correction(realPositionInRead, ErrorType::MULTIDEL, rwa.seq[realPositionInRead]));
}

void handleSoftClipping(ReadWithAlignments& rwa, size_t cigarCount, HandlingInfo& info) {
	for (size_t j = 0; j < cigarCount; ++j) {
		rwa.seq[info.positionInRead + info.hardClippedBases + j] = 'S';
	}
	info.positionInRead += cigarCount;
	info.softClippedBases += cigarCount;
}

void handleInsertion(ReadWithAlignments& rwa, size_t cigarCount, HandlingInfo& info) {
	for (size_t j = 0; j < cigarCount; ++j) {
		unsigned nucleotidePositionRead = info.positionInRead;
		char nucleotideInRead = rwa.seq[nucleotidePositionRead];

		if (nucleotideInRead == 'S') {
			throw std::runtime_error("this should not happen");
		}

		/*unsigned nucleotidePositionInReference = info.beginPos + nucleotidePositionRead - info.insertedBases
		 - info.softClippedBases + info.deletedBases;*/
		size_t realPositionInRead = info.positionInRead + info.hardClippedBases;

		if (info.revComp) {
			nucleotideInRead = util::reverseComplementChar(nucleotideInRead);
		}

		assert(realPositionInRead < readLength);

		char fromBase = nucleotideInRead;

		size_t correctionPosition = realPositionInRead;
		if (info.revComp) {
			correctionPosition = rwa.seq.size() - correctionPosition - 1;
		}

		info.corrections.push_back(Correction(correctionPosition, ErrorType::INSERTION, fromBase));
	}
	info.insertedBases += cigarCount;
}

void handleDeletion(ReadWithAlignments& rwa, size_t cigarCount, HandlingInfo& info, const std::string& genome) {
	size_t realPositionInRead = info.positionInRead + info.hardClippedBases - 1;
	assert(realPositionInRead < readLength);
	char fromBase = rwa.seq[realPositionInRead];
	if (fromBase == 'S') {
		throw std::runtime_error("this should not happen");
	}

	std::string toBases = "";
	size_t correctionPosition = realPositionInRead;
	if (info.revComp) {
		correctionPosition = rwa.seq.size() - correctionPosition - 1;
	}

	unsigned nucleotidePositionInReference = (info.positionInRead + info.hardClippedBases) + info.beginPos
			- info.softClippedBases;

	if (nucleotidePositionInReference >= genome.size()) {
		throw std::runtime_error("nucleotidePositionInReference >= length(genome)");
	}

	char nucleotideInReference;
	for (size_t j = 0; j < cigarCount; ++j) {
		info.deletedBases++;
		nucleotideInReference = genome[nucleotidePositionInReference + j];
		if (info.revComp) {
			nucleotideInReference = util::reverseComplementChar(nucleotideInReference);
		}
		toBases += nucleotideInReference;
	}
	assert(toBases.size() == cigar[i].count);
	if (cigarCount == 1) {
		if (nucleotideInReference == 'A') {
			info.corrections.push_back(Correction(correctionPosition, ErrorType::DEL_OF_A, fromBase));
		} else if (nucleotideInReference == 'C') {
			info.corrections.push_back(Correction(correctionPosition, ErrorType::DEL_OF_C, fromBase));
		} else if (nucleotideInReference == 'G') {
			info.corrections.push_back(Correction(correctionPosition, ErrorType::DEL_OF_G, fromBase));
		} else if (nucleotideInReference == 'T') {
			info.corrections.push_back(Correction(correctionPosition, ErrorType::DEL_OF_T, fromBase));
		} else if (nucleotideInReference != 'N') {
			throw std::runtime_error("weird base in reference genome");
		}
	} else {
		info.corrections.push_back(Correction(correctionPosition, ErrorType::MULTIDEL, fromBase));
	}
	info.positionInRead += cigarCount;
}

// Assumes that indels have already been detected
void handleSubstitutionErrors(ReadWithAlignments& rwa, HandlingInfo& info, const std::string& genomePath) {
	int genomeIdx = info.beginPos - 1;

	const std::string genome = io::readReferenceGenome(genomePath);

	for (size_t i = 0; i < rwa.seq.size(); ++i) {
		if (rwa.seq[i] == 'S') {
			continue;
		} else {
			genomeIdx++;
		}

		char baseInRead = rwa.seq[i];
		char baseInGenome = genome[genomeIdx];

		size_t correctionPos = i;

		if (info.revComp) {
			correctionPos = rwa.seq.size() - i - 1;
			baseInRead = util::reverseComplementChar(rwa.seq[correctionPos]);
		}

		if (baseInRead != baseInGenome) {
			std::string debugString = "";
			for (size_t t = 0; t < rwa.seq.size(); ++t) {
				debugString += genome[info.beginPos + t];
			}

			// TODO: Check me.
			if (baseInGenome == 'A') {
				if (!info.revComp) {
					info.corrections.push_back(Correction(correctionPos, ErrorType::SUB_OF_A, baseInRead));
				} else {
					info.corrections.push_back(Correction(correctionPos, ErrorType::SUB_OF_T, baseInRead));
				}
			} else if (baseInGenome == 'C') {
				if (!info.revComp) {
					info.corrections.push_back(Correction(correctionPos, ErrorType::SUB_OF_C, baseInRead));
				} else {
					info.corrections.push_back(Correction(correctionPos, ErrorType::SUB_OF_G, baseInRead));
				}
			} else if (baseInGenome == 'G') {
				if (!info.revComp) {
					info.corrections.push_back(Correction(correctionPos, ErrorType::SUB_OF_G, baseInRead));
				} else {
					info.corrections.push_back(Correction(correctionPos, ErrorType::SUB_OF_C, baseInRead));
				}
			} else if (baseInGenome == 'T') {
				if (!info.revComp) {
					info.corrections.push_back(Correction(correctionPos, ErrorType::SUB_OF_T, baseInRead));
				} else {
					info.corrections.push_back(Correction(correctionPos, ErrorType::SUB_OF_A, baseInRead));
				}
			} else if (baseInGenome == 'N') {
				throw std::runtime_error("base in genome is N");
			} else {
				std::cout << baseInGenome << "\n";
				throw std::runtime_error("Weird base in genome");
			}
		}
	}
	rwa.endPos = genomeIdx - 1;
}

std::vector<Correction> extractErrors(ReadWithAlignments& rwa, const std::string &genome) {
	if (hasFlagUnmapped(rwa.records[0])) {
		throw std::runtime_error("The read is unmapped");
	}
	// warn about searching errors in soft-clipped non-chimeric reads
	if (rwa.records.size() == 1) {
		for (size_t i = 0; i < length(rwa.records[0].cigar); ++i) {
			if (rwa.records[0].cigar[i].operation == 'S') {
				std::cout << "WARNING: Searching errors in soft-clipped non-chimeric read" << std::endl;
			}
		}
	}
	// warn about searching errors in reads with chimeric breaks
	if (rwa.records.size() > 1) {
		std::cout << "WARNING: chimeric break detected" << std::endl;
	}

	HandlingInfo info;
	for (seqan::BamAlignmentRecord record : rwa.records) {
		// Extract CIGAR string
		std::string cigarString = extractCigarString(record.cigar);
		info.positionInRead = 0;
		info.hardClippedBases = 0;
		info.softClippedBases = 0;
		info.insertedBases = 0;
		info.deletedBases = 0;
		info.readLength = rwa.seq.size();
		info.beginPos = record.beginPos;
		info.revComp = hasFlagRC(record);
		for (unsigned i = 0; i < length(record.cigar); ++i) {
			if (record.cigar[i].operation == 'H') { // hard clipping
				handleChimericBreak(rwa, record.cigar[i].count, info);
			} else if (record.cigar[i].operation == 'S') { // soft clipping
				handleSoftClipping(rwa, record.cigar[i].count, info);
			} else if (record.cigar[i].operation == 'M') { // match ... check for substitutions later
				info.positionInRead += record.cigar[i].count;
			} else if (record.cigar[i].operation == 'I') { // insertion
				handleInsertion(rwa, record.cigar[i].count, info);
			} else if (record.cigar[i].operation == 'D') { // deletion
				handleDeletion(rwa, record.cigar[i].count, info, genome);
			} else if (record.cigar[i].operation == 'P' || record.cigar[i].operation == 'N') {
				std::cout << "ERROR! There are letters I don't understand yet!" << record.cigar[i].operation
						<< std::endl;
			}
		}
		// now that indels have been fixed, fix the substitution errors.
		handleSubstitutionErrors(rwa, info, genome);
	}
	std::sort(info.corrections.begin(), info.corrections.end());
	return info.corrections;
}

template<typename T>
std::vector<std::vector<T> > createMatrix(size_t nrows, size_t ncols) {
	std::vector<std::vector<T> > res;
	res.resize(nrows);
	for (size_t i = 0; i < nrows; ++i) {
		res[i].resize(ncols);
	}
	return res;
}

std::vector<std::vector<int> > fillDPMatrix(const std::string& s1, const std::string& s2) {
	// create a DP table
	size_t nrows = s1.size() + 1;
	size_t ncols = s2.size() + 1;
	std::vector<std::vector<int> > matrix = createMatrix<int>(nrows, ncols);
	// initialize it
	for (size_t i = 0; i < nrows; ++i) {
		matrix[i][0] = i;
	}
	for (size_t j = 0; j < ncols; ++j) {
		matrix[0][j] = j;
	}
	// fill it
	for (size_t i = 1; i < nrows; ++i) {
		for (size_t j = 1; j < ncols; ++j) {
			int penalty = (s1[i] == s2[j]) ? 0 : 1;
			matrix[i][j] = std::min(matrix[i - 1][j - 1] + penalty,
					std::min(matrix[i - 1][j] + 1, matrix[i][j - 1] + 1));
		}
	}
	return matrix;
}

/*
 std::vector<AlignedCorrection> extractErrors(const std::string& originalRead, const std::string& correctedRead) {
 std::vector<std::vector<int> > dpMatrix = fillDPMatrix(originalRead, correctedRead);


 throw std::runtime_error("not implemented yet");
 }
 */

void handleUndetectedError(size_t posTruth, ErrorType typeTruth, EvaluationData& data, std::vector<bool>& fineBases,
		std::vector<bool>& fineGaps) {
	if (isBaseErrorType(typeTruth)) {
		data.update(typeTruth, ErrorType::CORRECT);
		fineBases[posTruth] = false;
	} else {
		data.update(typeTruth, ErrorType::NODEL);
		fineGaps[posTruth] = false;
	}
}

void handleMisdetectedError(size_t posPredicted, ErrorType typePredicted, EvaluationData& data,
		std::vector<bool>& fineBases, std::vector<bool>& fineGaps) {
	if (isBaseErrorType(typePredicted)) {
		data.update(ErrorType::CORRECT, typePredicted);
		fineBases[posPredicted] = false;
	} else {
		data.update(ErrorType::NODEL, typePredicted);
		fineGaps[posPredicted] = false;
	}
}

// requires the detected errors to be sorted by position in read
void updateEvaluationData(EvaluationData& data, const std::vector<Correction>& errorsTruth,
		const std::vector<Correction>& errorsPredicted, size_t readLength, const std::string& mappedSequence) {
	size_t truthIdx = 0;
	size_t predictedIdx = 0;

	std::vector<bool> fineBases(readLength, true);
	std::vector<bool> fineGaps(readLength, true);

	// Reißverschlussverfahren
	while (truthIdx < errorsTruth.size() && predictedIdx < errorsPredicted.size()) {
		size_t posTruth = errorsTruth[truthIdx].positionInRead;
		size_t posPredicted = errorsPredicted[predictedIdx].positionInRead;
		ErrorType typeTruth = errorsTruth[truthIdx].errorType;
		ErrorType typePredicted = errorsPredicted[predictedIdx].errorType;

		if (posTruth < posPredicted) { // undetected error
			if (mappedSequence[posTruth] != 'S') {
				handleUndetectedError(posTruth, typeTruth, data, fineBases, fineGaps);
			}
			truthIdx++;
		} else if (posTruth == posPredicted) { // errors at same position
			if (isBaseErrorType(typeTruth) && isBaseErrorType(typePredicted)) {
				if (mappedSequence[posTruth] != 'S') {
					data.update(typeTruth, typePredicted);
					fineBases[posTruth] = false;
				}
				truthIdx++;
				predictedIdx++;
			} else if (isGapErrorType(typeTruth) && isGapErrorType(typePredicted)) {
				if (mappedSequence[posTruth] != 'S') {
					data.update(typeTruth, typePredicted);
					fineGaps[posTruth] = false;
				}
				truthIdx++;
				predictedIdx++;
			} else {
				if (typeTruth < typePredicted) { // undetected error
					if (mappedSequence[posTruth] != 'S') {
						handleUndetectedError(posTruth, typeTruth, data, fineBases, fineGaps);
					}
					truthIdx++;
				} else { // misdetected error
					if (mappedSequence[posPredicted] != 'S') {
						handleMisdetectedError(posPredicted, typePredicted, data, fineBases, fineGaps);
					}
					predictedIdx++;
				}
			}
		} else { // misdetected error
			if (mappedSequence[posPredicted] != 'S') {
				handleMisdetectedError(posPredicted, typePredicted, data, fineBases, fineGaps);
			}
			predictedIdx++;
		}
	}
	while (truthIdx < errorsTruth.size()) { // further undetected errors
		size_t posTruth = errorsTruth[truthIdx].positionInRead;
		ErrorType typeTruth = errorsTruth[truthIdx].errorType;
		if (mappedSequence[posTruth] != 'S') {
			handleUndetectedError(posTruth, typeTruth, data, fineBases, fineGaps);
		}
		truthIdx++;
	}
	while (predictedIdx < errorsPredicted.size()) { // further misdetected errors
		size_t posPredicted = errorsPredicted[predictedIdx].positionInRead;
		ErrorType typePredicted = errorsPredicted[predictedIdx].errorType;
		if (mappedSequence[posPredicted] != 'S') {
			handleMisdetectedError(posPredicted, typePredicted, data, fineBases, fineGaps);
		}
		predictedIdx++;
	}

	// correctly unchanged bases/gaps
	for (size_t i = 0; i < readLength; ++i) {
		if (mappedSequence[i] == 'S') {
			continue;
		}
		if (fineBases[i]) {
			data.update(ErrorType::CORRECT, ErrorType::CORRECT);
		}
		if (fineGaps[i]) {
			data.update(ErrorType::NODEL, ErrorType::NODEL);
		}
	}
}

std::vector<Correction> convertToCorrections(const std::vector<std::pair<size_t, ErrorType> >& alignment,
		const std::string& originalRead) {
	std::vector<Correction> res;

	for (size_t i = 0; i < alignment.size(); ++i) {
		res.push_back(Correction(alignment[i].first, alignment[i].second, originalRead[alignment[i].first]));
	}

	return res;
}

EvaluationData evaluateCorrectionsByAlignment(const std::string& alignedReadsFilepath,
		const std::string& correctedReadsFilepath, const std::string& genomeFilepath) {
	EvaluationData data;
	//std::string genome = io::readReferenceGenome(genomeFilepath);
	BAMIterator it(alignedReadsFilepath);
	std::ifstream infile(correctedReadsFilepath);

	io::ReadInput input;
	input.openFile(correctedReadsFilepath);

	while (it.hasReadsLeft() && input.hasNext()) {
		io::Read correctedRead = input.readNext(true, true, true);
		ReadWithAlignments rwa = it.next();
		while (hasFlagUnmapped(rwa.records[0])) {
			rwa = it.next();
			correctedRead = input.readNext(true, true, true);
		}
		if (correctedRead.name.substr(0, correctedRead.name.find('/')) != rwa.name) {
			throw std::runtime_error("Something went wrong while evaluating the reads. Are they really sorted?\n Corrected read name: " + correctedRead.name + "\nOriginal read name: " + rwa.name);
		}
		std::vector<Correction> errorsTruth = extractErrors(rwa, genomeFilepath);
		std::vector<Correction> errorsPredicted = convertToCorrections(align(rwa.seq, correctedRead.seq), rwa.seq);
		updateEvaluationData(data, errorsTruth, errorsPredicted, correctedRead.seq.size(), rwa.seq);
	}
	return data;
}

} // end of namespace seq_correct::eval
} // end of namespace seq_correct
