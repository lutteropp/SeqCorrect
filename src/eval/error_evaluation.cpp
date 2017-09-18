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
#include "../util/util.hpp"
#include "../io/bam_iterator.hpp"
#include "../io/read_with_alignments.hpp"
#include "../io/sequence_io.hpp"
#include "alignment.hpp"
#include "../kmer/classification.hpp"
#include "../kmer/bloom_filter_classifier.hpp"
#include "error_evaluation.hpp"

namespace seq_correct {
namespace eval {

using namespace classification;

void printErrorEvaluationData(const eval::ErrorEvaluationData& evalData);

// Compare the corrected reads with the aligned original reads
/*EvaluationData evaluateCorrections(const std::string& originalReadsFilepath, const std::string& correctedReadsFilepath,
 const std::string& genomeFilepath);*/

ErrorEvaluationData evaluateCorrectionsByAlignment(const std::string& alignmentFilepath,
		const std::string& correctedReadsFilepath, const std::string& genomeFilepath, util::GenomeType genomeType);

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

void handleChimericBreak(std::string& seq, size_t cigarCount, HandlingInfo& info) {
	/*unsigned nucleotidePositionInReference = info.beginPos + info.positionInRead - info.insertedBases
	 - info.softClippedBases + info.deletedBases;*/
	unsigned realPositionInRead = info.positionInRead + info.hardClippedBases; // the position of the first base of the chimeric break

	assert(realPositionInRead < seq.size());

	if (realPositionInRead == 0) { // in this case, take the position of the last base of the chimeric break
		realPositionInRead += cigarCount;
	}

	size_t correctionPosition = realPositionInRead;
	if (info.revComp) {
		correctionPosition = seq.size() - realPositionInRead - 1;
	}

	// What if multiple chimeric breaks happen in a read? ... It should be fine, too.
	//  As long as there aren't multiple hard clipped regions in a single record.
	info.hardClippedBases += cigarCount;

	// because we treat chimeric breaks as deletion of multiple bases
	info.corrections.push_back(Correction(correctionPosition, ErrorType::MULTIDEL, seq[realPositionInRead]));
}

void handleSoftClipping(std::string& seq, size_t cigarCount, HandlingInfo& info) {
	if (info.revComp) {
		seq = reverseComplementString(seq);
	}

	for (size_t j = 0; j < cigarCount; ++j) {
		seq[info.positionInRead + info.hardClippedBases + j] = 'S';
	}
	info.positionInRead += cigarCount;
	info.softClippedBases += cigarCount;

	if (info.revComp) {
		seq = reverseComplementString(seq);
	}
}

void handleInsertion(std::string& seq, size_t cigarCount, HandlingInfo& info) {
	if (info.revComp) {
		seq = reverseComplementString(seq);
	}

	for (size_t j = 0; j < cigarCount; ++j) {
		unsigned nucleotidePositionRead = info.positionInRead;
		char nucleotideInRead = seq[nucleotidePositionRead];

		if (nucleotideInRead == 'S') {
			std::cout << "positionInRead: " << info.positionInRead << "\n";
			std::cout << "seq: " << seq << "\n";
			std::cout << "insertion error detected, but the base is an S\n";
			throw std::runtime_error("this should not happen");
		}

		size_t realPositionInRead = info.positionInRead + info.hardClippedBases;
		assert(realPositionInRead < seq.size());

		if (realPositionInRead > seq.size()) {
			throw std::runtime_error(
					"realPositionInRead = " + std::to_string(realPositionInRead) + ", but seq.size() = "
							+ std::to_string(seq.size()));
		}

		size_t correctionPosition = realPositionInRead;
		char fromBase = nucleotideInRead;
		if (info.revComp) {
			correctionPosition = seq.size() - realPositionInRead - 1;
			if (correctionPosition > seq.size()) {
				throw std::runtime_error(
						"correctionPosition = " + std::to_string(correctionPosition) + ", but seq.size() = "
								+ std::to_string(seq.size()));
			}
			fromBase = reverseComplementChar(fromBase);
		}
		info.corrections.push_back(Correction(correctionPosition, ErrorType::INSERTION, fromBase));

		seq[realPositionInRead] = tolower(seq[realPositionInRead]); // mark the insertion
		//rwa.seq = rwa.seq.substr(0, correctionPosition) + rwa.seq.substr(correctionPosition + 1, std::string::npos);
	}
	info.insertedBases += cigarCount;

	if (info.revComp) {
		seq = reverseComplementString(seq);
	}
}

void handleDeletion(std::string& seq, size_t cigarCount, HandlingInfo& info, const std::string& genome) {
	if (info.revComp) {
		seq = reverseComplementString(seq);
	}

	size_t realPositionInRead = info.positionInRead + info.hardClippedBases - 1;
	assert(realPositionInRead < seq.size());
	char fromBase = seq[realPositionInRead];
	/*if (fromBase == 'S') {
	 std::cout << "deletion error detected, but the base is an S\n";
	 throw std::runtime_error("this should not happen");
	 }*/

	unsigned nucleotidePositionInReference = (info.positionInRead + info.hardClippedBases) + info.beginPos
			- info.softClippedBases;

	char nucleotideInReference;
	for (size_t j = 0; j < cigarCount; ++j) {
		info.deletedBases++;
		nucleotideInReference = genome[nucleotidePositionInReference + j];
	}

	size_t correctionPosition = realPositionInRead;
	if (info.revComp) {
		correctionPosition = seq.size() - realPositionInRead - 1;
	}
	if (info.revComp) {
		fromBase = reverseComplementChar(fromBase);
	}

	if (cigarCount == 1) {
		if (nucleotideInReference == 'A') {
			if (!info.revComp) {
				info.corrections.push_back(Correction(correctionPosition, ErrorType::DEL_OF_A, fromBase));
			} else {
				info.corrections.push_back(Correction(correctionPosition, ErrorType::DEL_OF_T, fromBase));
			}
		} else if (nucleotideInReference == 'C') {
			if (!info.revComp) {
				info.corrections.push_back(Correction(correctionPosition, ErrorType::DEL_OF_C, fromBase));
			} else {
				info.corrections.push_back(Correction(correctionPosition, ErrorType::DEL_OF_G, fromBase));
			}
		} else if (nucleotideInReference == 'G') {
			if (!info.revComp) {
				info.corrections.push_back(Correction(correctionPosition, ErrorType::DEL_OF_G, fromBase));
			} else {
				info.corrections.push_back(Correction(correctionPosition, ErrorType::DEL_OF_C, fromBase));
			}
		} else if (nucleotideInReference == 'T') {
			if (!info.revComp) {
				info.corrections.push_back(Correction(correctionPosition, ErrorType::DEL_OF_T, fromBase));
			} else {
				info.corrections.push_back(Correction(correctionPosition, ErrorType::DEL_OF_A, fromBase));
			}
		} else if (nucleotideInReference != 'N') {
			throw std::runtime_error("weird base in reference genome");
		}
	} else {
		info.corrections.push_back(Correction(correctionPosition, ErrorType::MULTIDEL, fromBase));
	}

	std::string deletedBases = "";
	for (size_t i = 0; i < cigarCount; ++i) {
		deletedBases += "D";
	}
	seq = seq.substr(0, realPositionInRead + 1) + deletedBases + seq.substr(realPositionInRead + 1, std::string::npos);

	info.positionInRead += cigarCount;
	info.deletedBases += cigarCount;

	if (info.revComp) {
		seq = reverseComplementString(seq);
	}
}

// Assumes that indels have already been detected
void handleSubstitutionErrors(std::string& seq, HandlingInfo& info, const std::string& genome, GenomeType genomeType) {
	if (info.revComp) {
		seq = reverseComplementString(seq);
	}
	int genomeIdx = -1;

	std::string genomeDebugString = genome.substr(info.beginPos, seq.size());

	for (size_t i = 0; i < seq.size(); ++i) {
		if (seq[i] == 'S' || islower(seq[i])) {
			continue;
		} else {
			genomeIdx++;
		}

		char baseInRead = seq[i];
		if (baseInRead == 'D') { // marked position of a deletion error
			continue;
		}

		char baseInGenome;
		if (genomeType == GenomeType::LINEAR) {
			baseInGenome = genome[info.beginPos + genomeIdx];
		} else {
			baseInGenome = genome[(info.beginPos + genomeIdx) % genome.size()];
		}

		size_t correctionPos = i;

		if (baseInRead != baseInGenome) {
			if (info.revComp) {
				correctionPos = seq.size() - i - 1;
				baseInRead = reverseComplementChar(baseInRead);
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

	if (info.revComp) {
		seq = reverseComplementString(seq);
	}
}

void fixReadBackToNormal(ReadWithAlignments& rwa) {
	std::string seq = "";
	for (size_t i = 0; i < rwa.seq.size(); ++i) {
		if (rwa.seq[i] != 'D') {
			seq += toupper(rwa.seq[i]);
		}
	}
	rwa.seq = seq;
}

std::vector<Correction> extractErrors(const ReadWithAlignments& rwa, const std::string &genome, GenomeType genomeType) {
	if (hasFlagUnmapped(rwa.records[0])) {
		throw std::runtime_error("The read is unmapped");
	}
	// warn about searching errors in soft-clipped non-chimeric reads
	/*if (rwa.records.size() == 1) {
	 for (size_t i = 0; i < length(rwa.records[0].cigar); ++i) {
	 if (rwa.records[0].cigar[i].operation == 'S') {
	 std::cout << "WARNING: Searching errors in soft-clipped non-chimeric read" << std::endl;
	 }
	 }
	 }*/
	// warn about searching errors in reads with chimeric breaks
	if (rwa.records.size() > 1) {
		std::cout << "WARNING: chimeric break detected" << std::endl;
	}

	HandlingInfo info;
	std::string seq = rwa.seq;

	if (hasFlagRC(rwa.records[0])) {
		info.revComp = true;
	}

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
		for (unsigned i = 0; i < length(record.cigar); ++i) {
			if (record.cigar[i].operation == 'H') { // hard clipping
				handleChimericBreak(seq, record.cigar[i].count, info);
			} else if (record.cigar[i].operation == 'S') { // soft clipping
				handleSoftClipping(seq, record.cigar[i].count, info);

				if (seq[info.positionInRead] == 'S') {
					std::cout << "cigar: " << cigarString << "\n";
					std::cout << "original: " << rwa.seq << "\n";
					std::cout << "current seq: " << seq << "\n";
					std::cout << "current posInRead: " << info.positionInRead << "\n";
				}

			} else if (record.cigar[i].operation == 'M') { // match ... check for substitutions later
				info.positionInRead += record.cigar[i].count;
			} else if (record.cigar[i].operation == 'I') { // insertion

				if (seq == "CCAGATCAGAGTTTTGTTTAGGAGCTGGGTCCTCCCTGATGSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS") {
					std::cout << "I'm here for debugging.\n";
					std::cout << "cigar: " << cigarString << "\n";
					std::cout << "orig: " << rwa.seq << "\n";
				}

				handleInsertion(seq, record.cigar[i].count, info);
			} else if (record.cigar[i].operation == 'D') { // deletion
				handleDeletion(seq, record.cigar[i].count, info, genome);
			} else if (record.cigar[i].operation == 'P' || record.cigar[i].operation == 'N') {
				std::cout << "ERROR! There are letters I don't understand yet!" << record.cigar[i].operation
						<< std::endl;
			}
		}
		// now that indels have been fixed, fix the substitution errors.
		handleSubstitutionErrors(seq, info, genome, genomeType);
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

void handleUndetectedError(size_t posTruth, ErrorType typeTruth, ErrorEvaluationData& data,
		std::vector<bool>& fineBases, std::vector<bool>& fineGaps) {
	if (isBaseErrorType(typeTruth)) {
		data.update(typeTruth, ErrorType::CORRECT);
		fineBases[posTruth] = false;
	} else {
		data.update(typeTruth, ErrorType::NODEL);
		fineGaps[posTruth] = false;
	}
}

void handleMisdetectedError(size_t posPredicted, ErrorType typePredicted, ErrorEvaluationData& data,
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
void updateEvaluationData(ErrorEvaluationData& data, const std::vector<Correction>& errorsTruth,
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

// do some fancy read-name-parsing to comply with 'samtools sort -n'
static int strnum_cmp(const std::string& s1, const std::string& s2) {
	// parse the strings into digits and non-digits
	std::vector<std::string> s1Parsed;
	std::vector<std::string> s2Parsed;
	bool s1Digit = isdigit(s1[0]);
	bool s2Digit = isdigit(s2[0]);
	bool startText1 = !s1Digit;
	bool startText2 = !s2Digit;
	s1Parsed.push_back(std::string(""));
	s1Parsed[0] += s1[0];
	s2Parsed.push_back(std::string(""));
	s2Parsed[0] += s2[0];
	for (size_t i = 1; i < s1.size(); ++i) {
		if (isdigit(s1[i]) == s1Digit) {
			s1Parsed[s1Parsed.size() - 1] += s1[i];
		} else {
			s1Digit = !s1Digit;
			s1Parsed.push_back(std::string(""));
			s1Parsed[s1Parsed.size() - 1] += s1[i];
		}
	}
	for (size_t i = 1; i < s2.size(); ++i) {
		if (isdigit(s2[i]) == s2Digit) {
			s2Parsed[s2Parsed.size() - 1] += s2[i];
		} else {
			s2Digit = !s2Digit;
			s2Parsed.push_back(std::string(""));
			s2Parsed[s2Parsed.size() - 1] += s2[i];
		}
	}

	if (startText1 != startText2) {
		return (s1Parsed[0] < s2Parsed[0]);
	}
	for (size_t i = 0; i < std::min(s1Parsed.size(), s2Parsed.size()); ++i) {
		if (s1Parsed[i] != s2Parsed[i]) {
			if (isdigit(s1Parsed[i][0]) && isdigit(s2Parsed[i][0])) {
				// compare numbers
				return (stoi(s1Parsed[i]) < stoi(s2Parsed[i]));
			} else {
				return (s1Parsed[i] < s2Parsed[i]);
			}
		}
	}
	return (s1Parsed.size() < s2Parsed.size());
}

// TODO: Maybe replace this by external Merge sort if the read data set is too large
void sortCorrectedReads(const std::string& correctedReadsFilepath, const std::string& sortedFilepath) {
	std::ifstream corrTest(sortedFilepath);
	if (corrTest.good()) {
		std::cout << "Sorted reads file already exists. Doing nothing." << "\n";
		return;
	}

	io::ReadOutput output;
	output.createFile(sortedFilepath);

	std::vector<io::Read> allReads;
	io::ReadInput input(correctedReadsFilepath);
	while (input.hasNext()) {
		allReads.push_back(input.readNext(true, false, true));
	}
	std::cout << allReads.size() << "\n";
	std::sort(allReads.begin(), allReads.end(),
			[](const io::Read& left, const io::Read& right) {return strnum_cmp(left.name, right.name);});
	for (size_t i = 0; i < allReads.size(); ++i) {
		output.write(allReads[i]);
	}
	output.close();
}

ErrorEvaluationData evaluateCorrectionsByAlignment(const std::string& alignedReadsFilepath,
		const std::string& correctedReadsFilepath, const std::string& genomeFilepath, util::GenomeType genomeType) {
	ErrorEvaluationData data;
	std::string genome = io::readReferenceGenome(genomeFilepath);
	BAMIterator it(alignedReadsFilepath);

	std::ifstream test(correctedReadsFilepath);
	if (!test.good()) {
		throw std::runtime_error("The specified file does not exist: " + correctedReadsFilepath);
	}
	test.close();

	std::string sortedReadsFilepath = correctedReadsFilepath + ".sorted";
	sortCorrectedReads(correctedReadsFilepath, sortedReadsFilepath);

	io::ReadInput input(sortedReadsFilepath);

#pragma omp parallel
	{
#pragma omp single
		{
			while (it.hasReadsLeft() && input.hasNext()) {
				io::Read correctedRead = input.readNext(true, false, true);
				//std::cout << "correctedRead = " << correctedRead.name << "\n";

				ReadWithAlignments rwa = it.next();
				//std::cout << "rwa = " << rwa.name << "\n";
				while (hasFlagUnmapped(rwa.records[0]) && input.hasNext() && it.hasReadsLeft()) {
					std::cout << "skipping " << rwa.name << " because it's unmapped\n";
					rwa = it.next();
					//std::cout << "rwa = " << rwa.name << "\n";
					std::cout << "skipping " << correctedRead.name << "\n";
					correctedRead = input.readNext(true, false, true);
					//std::cout << "correctedRead = " << correctedRead.name << "\n";
				}

				while (rwa.records.size() > 1 && input.hasNext() && it.hasReadsLeft()) {
					std::cout << "skipping " << rwa.name << " because it's chimeric\n";
					rwa = it.next();
					//std::cout << "rwa = " << rwa.name << "\n";
					std::cout << "skipping " << correctedRead.name << "\n";
					correctedRead = input.readNext(true, false, true);
					//std::cout << "correctedRead = " << correctedRead.name << "\n";
				}

				while (correctedRead.name.substr(0, correctedRead.name.find(' ')).substr(0,
						correctedRead.name.find('/')) != rwa.name && input.hasNext()) {
					correctedRead = input.readNext(true, false, true);
					//std::cout << "correctedRead = " << correctedRead.name << "\n";
					/*throw std::runtime_error(
					 "Something went wrong while evaluating the reads. Are they really sorted?\n Corrected read name: "
					 + correctedRead.name + "\nOriginal read name: " + rwa.name);*/
				}

				if (hasFlagUnmapped(rwa.records[0]) || rwa.records.size() > 1
						|| correctedRead.name.substr(0, correctedRead.name.find(' ')).substr(0,
								correctedRead.name.find('/')) != rwa.name) {
					break;
				}

#pragma omp task shared(genome, data), firstprivate(rwa, correctedRead)
				{
					std::vector<Correction> errorsTruth = extractErrors(rwa, genome, genomeType);
					std::vector<Correction> errorsPredicted;
					errorsPredicted = convertToCorrections(align(rwa.seq, correctedRead.seq), rwa.seq);

					if (errorsTruth.size() > 0 && errorsPredicted.size() > 0) {
						if (errorsTruth[0].errorType != errorsPredicted[0].errorType) {

							std::cout << "Original: " << rwa.seq << "\n";
							std::cout << "Corrected : " << correctedRead.seq << "\n";
							std::string genomeArea = genome.substr(rwa.beginPos, rwa.endPos - rwa.beginPos + 1);
							if (hasFlagRC(rwa.records[0])) {
								genomeArea = reverseComplementString(genomeArea);
							}
							std::cout << "genomeArea: " << genomeArea << "\n";

							std::cout << "errorsTruth:\n";
							for (size_t i = 0; i < errorsTruth.size(); ++i) {
								std::cout << "  " << errorTypeToString(errorsTruth[i].errorType) << " at "
										<< errorsTruth[i].positionInRead << "\n";
							}
							std::cout << "errorsPredicted:\n";
							for (size_t i = 0; i < errorsPredicted.size(); ++i) {
								std::cout << "  " << errorTypeToString(errorsPredicted[i].errorType) << " at "
										<< errorsPredicted[i].positionInRead << "\n";
							}
							if (hasFlagRC(rwa.records[0])) {
								std::cout << "The read is reverse-complemented.\n";
							}

							/*errorsTruth = extractErrors(rwa, genome, genomeType);
							 if (hasFlagRC(rwa.records[0])) {
							 errorsPredicted = convertToCorrections(align(rwa.seq, correctedRead.seq),
							 util::reverseComplementString(rwa.seq));
							 } else {
							 errorsPredicted = convertToCorrections(align(rwa.seq, correctedRead.seq), rwa.seq);
							 }*/

						}
					}

					updateEvaluationData(data, errorsTruth, errorsPredicted, correctedRead.seq.size(), rwa.seq);
				}
			}
		}
#pragma omp taskwait
	}
	return data;
}

void printErrorEvaluationData(const eval::ErrorEvaluationData& evalData) {
	for (ErrorType type : AllErrorTypeIterator()) {
		std::cout << errorTypeToString(type) << ":\n";
		std::cout << "  TP:          " << eval::truePositives(type, evalData) << "\n";
		std::cout << "  TN:          " << eval::trueNegatives(type, evalData) << "\n";
		std::cout << "  FP:          " << eval::falsePositives(type, evalData) << "\n";
		std::cout << "  FN:          " << eval::falseNegatives(type, evalData) << "\n";
		std::cout << "  Accuracy:    " << eval::computeAccuracy(type, evalData) << "\n";
		std::cout << "  Precision:   " << eval::computePrecision(type, evalData) << "\n";
		std::cout << "  Recall:      " << eval::computeRecall(type, evalData) << "\n";
		std::cout << "  Specificity: " << eval::computeSpecificity(type, evalData) << "\n";
		std::cout << "  Sensitivity: " << eval::computeSensitivity(type, evalData) << "\n";
		std::cout << "  Gain:        " << eval::computeGain(type, evalData) << "\n";
		std::cout << "  F1-score:    " << eval::computeF1Score(type, evalData) << "\n";
	}
	std::cout << "\n";
	std::cout << "Unweighted Average Base F1-score: " << eval::computeUnweightedAverageBaseF1Score(evalData) << "\n";
	std::cout << "Base NMI-score:                   " << eval::computeBaseNMIScore(evalData) << "\n";
	std::cout << "Unweighted Average Gap F1-score:  " << eval::computeUnweightedAverageGapF1Score(evalData) << "\n";
	std::cout << "Gap NMI-score:                    " << eval::computeGapNMIScore(evalData) << "\n";

	std::cout << "\n Base Confusion Matrix:\n";
	for (ErrorType e1 : BaseTypeIterator()) {
		for (ErrorType e2 : BaseTypeIterator()) {
			double percentage = (double) evalData.getEntry(e1, e2) / evalData.sumTruth(e1);
			std::cout << "[" << util::errorTypeToString(e1) << "][" << util::errorTypeToString(e2) << "]: "
					<< evalData.getEntry(e1, e2) << " = " << percentage * 100 << "%" << "\n";
		}
	}

	std::cout << "\n Gap Confusion Matrix:\n";
	for (ErrorType e1 : GapTypeIterator()) {
		for (ErrorType e2 : GapTypeIterator()) {
			double percentage = (double) evalData.getEntry(e1, e2) / evalData.sumTruth(e1);
			std::cout << "[" << util::errorTypeToString(e1) << "][" << util::errorTypeToString(e2) << "]: "
					<< evalData.getEntry(e1, e2) << " = " << percentage * 100 << "%" << "\n";
		}
	}
}

// TODO: Maybe also provide the option to align the reads to the genome on-the-fly in this code, instead of calling another program?
void eval_corrections(size_t k, GenomeType genomeType, const std::string& pathToOriginalReads,
		const std::string& pathToCorrectedReads, const std::string& pathToGenome, const std::string& outputPath) {
	std::string alignmentPath = pathToOriginalReads.substr(0, pathToOriginalReads.find_last_of('.')) + ".bam";
	std::ifstream alignmentFile(alignmentPath);
	if (!alignmentFile.good()) {
		std::cout << "Alignment file " << alignmentPath << " not found. Building alignment now...\n";
		std::string indexGenomeCall = "bwa index " + pathToGenome;
		std::cout << "Calling: " << indexGenomeCall << "\n";
		int status = std::system(indexGenomeCall.c_str());
		if (!WIFEXITED(status)) {
			throw std::runtime_error("Something went wrong while indexing genome!");
		}
		std::string alignReadsCall = "bwa mem -L 999999999 " + pathToGenome + " " + pathToOriginalReads
				+ " > myreads_aln.sam";
		std::cout << "Calling: " << alignReadsCall << "\n";
		status = std::system(alignReadsCall.c_str());
		if (!WIFEXITED(status)) {
			throw std::runtime_error("Something went wrong while aligning reads!");
		}
		std::string removeUnmappedAndCompressCall = "samtools view -bS -F 4 myreads_aln.sam > myreads_aln.bam";
		std::cout << "Calling: " << removeUnmappedAndCompressCall << "\n";
		status = std::system(removeUnmappedAndCompressCall.c_str());
		if (!WIFEXITED(status)) {
			throw std::runtime_error("Something went wrong while removing unmapped reads and compressing the file!");
		}
		std::string sortCall = "samtools sort -n -o " + alignmentPath + " myreads_aln.bam";
		std::cout << "Calling: " << sortCall << "\n";
		status = std::system(sortCall.c_str());
		if (!WIFEXITED(status)) {
			throw std::runtime_error("Something went wrong while sorting the file!");
		}
		std::string removeCall1 = "rm myreads_aln.sam";
		status = std::system(removeCall1.c_str());
		if (!WIFEXITED(status)) {
			throw std::runtime_error("Something went wrong while removing temporary file myreads_aln.sam!");
		}
		std::string removeCall2 = "rm myreads_aln.bam";
		status = std::system(removeCall2.c_str());
		if (!WIFEXITED(status)) {
			throw std::runtime_error("Something went wrong while removing temporary file myreads_aln.bam!");
		}
	}
	alignmentFile.close();

	eval::ErrorEvaluationData res = eval::evaluateCorrectionsByAlignment(alignmentPath, pathToCorrectedReads,
			pathToGenome, genomeType);
	printErrorEvaluationData(res);
}

} // end of namespace seq_correct::eval
} // end of namespace seq_correct
