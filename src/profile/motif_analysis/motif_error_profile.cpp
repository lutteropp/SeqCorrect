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

#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

#include "motif_error_profile.hpp"


namespace seq_correct {
namespace profile {

MotifErrorProfile::MotifErrorProfile(Matcher &kmerCounter) :
		counter(kmerCounter) {
	finalized = false;
	motifTree = MotifTree(MAX_MOTIF_SIZE);
}

void MotifErrorProfile::updateMotifData(int positionInSequence, const std::string &sequence, const ErrorType &type) {
	for (int window = 0; window < MAX_MOTIF_SIZE; ++window) {
		for (int startPos = std::max(0, positionInSequence - window); startPos <= positionInSequence; ++startPos) {
			size_t endPos = startPos + window;
			if (endPos >= sequence.size()) {
				break;
			}
			int positionInMotif = positionInSequence - startPos;
			std::string motifString = sequence.substr(startPos, endPos - startPos + 1);
			motifTree[motifString].entries[positionInMotif].numErrors[type]++;
		}
	}
}

void MotifErrorProfile::check(const std::vector<Correction>& corrections, const Read &originalRead) {
	finalized = false;
	for (Correction corr : corrections) {
		updateMotifData(corr.positionInRead, originalRead.seq, corr.errorType);
		assert(corr.type != ErrorType::CORRECT && corr.type != ErrorType::NODEL);
	}
}

std::unordered_map<ErrorType, double> MotifErrorProfile::getErrorProbabilitiesFinalized(const Read &read,
		size_t positionInRead) {
	return getErrorProbabilitiesFinalized(read.seq, positionInRead);
}

std::unordered_map<ErrorType, double> MotifErrorProfile::getErrorProbabilitiesFinalized(const std::string &kmer,
		size_t positionInKmer) {
	assert(finalized);
	std::unordered_map<ErrorType, double> probs;
	for (ErrorType type : AllErrorTypeIterator()) {
		probs[type] = findMostSignificantZScore(type, kmer, positionInKmer);
	}

	if (kmer[positionInKmer] == 'A') {
		probs[ErrorType::SUB_OF_A] = std::numeric_limits<double>::lowest();
	} else if (kmer[positionInKmer] == 'C') {
		probs[ErrorType::SUB_OF_C] = std::numeric_limits<double>::lowest();
	} else if (kmer[positionInKmer] == 'G') {
		probs[ErrorType::SUB_OF_G] = std::numeric_limits<double>::lowest();
	} else if (kmer[positionInKmer] == 'T') {
		probs[ErrorType::SUB_OF_T] = std::numeric_limits<double>::lowest();
	}
	probs[ErrorType::CORRECT] = 0;
	probs[ErrorType::NODEL] = 0;

	return probs;
}

std::unordered_map<ErrorType, double> MotifErrorProfile::getErrorProbabilities(const Read &read,
		size_t positionInRead) {
	if (read.seq[positionInRead] == '_') {
		throw std::runtime_error("Invalid k-mer!");
	}
	if (finalized) {
		return getErrorProbabilitiesFinalized(read, positionInRead);
	} else {
		// TODO: Is this necessary?
		finalize();
		return getErrorProbabilitiesFinalized(read, positionInRead);
	}
}

std::unordered_map<ErrorType, double> MotifErrorProfile::getKmerErrorProbabilities(const std::string &kmer,
		size_t positionInKmer) {
	if (kmer.find("_") != std::string::npos) {
		throw std::runtime_error("Invalid k-mer!");
	}
	if (finalized) {
		return getErrorProbabilitiesFinalized(kmer, positionInKmer);
	} else {
		// TODO: Is this necessary?
		finalize();
		return getErrorProbabilitiesFinalized(kmer, positionInKmer);
	}
}

void MotifErrorProfile::storeErrorProfile(const std::string &filepath) {
	std::ofstream outfile(filepath, std::ios::binary);
	cereal::BinaryOutputArchive oarchive(outfile);
	oarchive(*this);
}

void MotifErrorProfile::loadErrorProfile(const std::string &filepath, Matcher &counter) {
	std::ifstream infile(filepath, std::ios::binary);
	if (!infile.good()) {
		throw std::runtime_error("The file " + filepath + " does not exist!");
	}
	cereal::BinaryInputArchive iarchive(infile);
	MotifErrorProfile mep(counter);
	iarchive(mep);
	motifTree = mep.motifTree;
	finalized = mep.finalized;
}

void MotifErrorProfile::plotErrorProfile() {

	for (ErrorType type : AllErrorTypeIterator()) {
		for (size_t i = 0; i < motifTree.size(); ++i) {
			std::string motifString = motifTree.indexToMotif(i);
			// TODO: Decide on whether the whole motif error profile should be printed or just the significant parts of it
			for (size_t j = 0; j < motifTree[i].entries.size(); ++j) {
				if ((motifTree[i].entries[j].zScore[type] >= 1) || (motifTree[i].entries[j].zScore[type] <= -1)) {
					std::cout << motifTree[i].entries[j].zScore[type] << "; " << "Motif " << motifString
							<< " with position: " << j << " for type " << errorTypeToString(type) << "\n";
				}
			}
		}
	}
}

void MotifErrorProfile::reset() {
	finalized = false;
	motifTree.reset();
}

void MotifErrorProfile::finalize() {
	if (finalized) {
		return;
	}
	computeZScores();
	finalized = true;
}

void MotifErrorProfile::computeZScores() {
	std::cout << "Computing Z Scores...\n";
	// The Z score only makes sense for a motif of at least length 3
	for (size_t l = 3; l <= MAX_MOTIF_SIZE; ++l) {
		size_t lengthOffset = 0;
		for (size_t j = 1; j < l; ++j) {
			lengthOffset += std::pow(5, j);
		}
		for (size_t i = 0; i < std::pow(5, l); ++i) {
			size_t motifIndex = lengthOffset + i;
			std::string motifString = motifTree.indexToMotif(motifIndex);
			std::string motifLeftString(motifString);
			motifLeftString.pop_back(); // w_1 ... w_{m-1}
			std::string motifRightString(motifString);
			motifRightString.erase(0, 1); // w_2 ... w_m
			std::string motifInnerString(motifLeftString);
			motifInnerString.erase(0, 1); // w_2... w_{m-1}
			size_t motifLeftIndex = motifTree.motifToIndex(motifLeftString);
			size_t motifRightIndex = motifTree.motifToIndex(motifRightString);
			size_t motifInnerIndex = motifTree.motifToIndex(motifInnerString);

			for (size_t pos = 0; pos < l; pos++) {
				// If pos is at an end of a motif, we need some k-mer counts. The main formula doesn't change. (Thanks to Nick Goldman for explaining this to me)
				for (ErrorType type : AllErrorTypeIterator()) {
					double countMotifInner;
					if ((pos == 0) || (pos == l - 1)) {
						countMotifInner = counter.countKmer(motifInnerString);
					} else {
						countMotifInner = motifTree[motifInnerIndex].entries[pos - 1].numErrors[type];
					}

					if (countMotifInner == 0) {
						continue;
					}

					double countMotif = motifTree[motifIndex].entries[pos].numErrors[type];
					double countMotifRight;
					if (pos == 0) {
						countMotifRight = counter.countKmer(motifRightString);
					} else {
						countMotifRight = motifTree[motifRightIndex].entries[pos - 1].numErrors[type];
					}

					double countMotifLeft;
					if (pos == l - 1) {
						countMotifLeft = counter.countKmer(motifLeftString);
					} else {
						countMotifLeft = motifTree[motifLeftIndex].entries[pos].numErrors[type];
					}

					double expectedCount = (countMotifLeft * countMotifRight) / countMotifInner;

					double variance = (countMotifInner - countMotifLeft) * (countMotifInner - countMotifRight)
							/ (countMotifInner * countMotifInner);
					variance *= expectedCount;

					if (variance != 0) {
						motifTree[motifIndex].entries[pos].zScore[type] = (countMotif - expectedCount) / sqrt(variance);
						assert(motifTree[motifIndex].entries[pos].zScore[type] > -999999999);
						assert(motifTree[motifIndex].entries[pos].zScore[type] < 999999999);
					}
				}
			}
		}
	}
	std::cout << "Finished computing Z Scores\n";
}

double MotifErrorProfile::findMostSignificantZScore(const ErrorType &type, const std::string &sequence,
		int posInSequence) {
	if (type == ErrorType::CORRECT || type == ErrorType::NODEL) {
		throw std::runtime_error("type is CORRECT or NODEL!");
	}
	if (!finalized) {
		finalize();
	}
	double minZScore = std::numeric_limits<double>::max();
	double maxZScore = std::numeric_limits<double>::lowest();
	for (int l = 3; l <= MAX_MOTIF_SIZE; ++l) {
		// build all possible motifs of size l from longest
		for (int start = std::max(0, posInSequence - l + 1); start <= posInSequence; start++) {
			size_t end = start + l - 1;
			if (end >= sequence.size()) {
				break;
			}
			std::string actMotif = sequence.substr(start, l);
			int posInActMotif = posInSequence - start;
			double actZScore = motifTree[actMotif].entries[posInActMotif].zScore[type];
			minZScore = std::min(minZScore, actZScore);
			maxZScore = std::max(maxZScore, actZScore);
		}
	}
	return (abs(minZScore) > abs(maxZScore)) ? minZScore : maxZScore;
}

} // end of namespace seq_correct::profile
} // end of namespace seq_correct
