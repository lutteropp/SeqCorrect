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
#include <functional>
#include <algorithm>
#include <vector>
#include "error_correction.hpp"

namespace seq_correct {
namespace correction {

using namespace classification;
using namespace std::placeholders;

Read correctRead_none(const Read& read, CorrectionParameters& params);
Read correctRead_simple_kmer(const Read& read, CorrectionParameters& params);
Read correctRead_adaptive_kmer(const Read& read, CorrectionParameters& params);
Read correctRead_full_msa(const Read& read, CorrectionParameters& params);
Read correctRead_partial_msa(const Read& read, CorrectionParameters& params);
Read correctRead_suffix_tree(const Read& read, CorrectionParameters& params);

std::pair<size_t, size_t> affectedReadArea(size_t posInRead, size_t readLength, size_t kmerSize) {
	size_t first = 0;
	if (posInRead + 1 > kmerSize) {
		first = posInRead - kmerSize + 1;
	}
	size_t last = std::min(readLength - 1, posInRead + kmerSize - 1);
	return std::make_pair(first, last);
}

std::string findReplacement(const std::string& kmer, ErrorType errorType, size_t posInKmer) {
	std::string replacement;
	if (errorType == ErrorType::SUB_OF_A) {
		replacement = "A";
	} else if (errorType == ErrorType::SUB_OF_C) {
		replacement = "C";
	} else if (errorType == ErrorType::SUB_OF_G) {
		replacement = "G";
	} else if (errorType == ErrorType::SUB_OF_T) {
		replacement = "T";
	} else if (errorType == ErrorType::INSERTION) {
		replacement = ""; // TODO: This will lead to k-mers of even size. Is this a problem?
	} else {
		replacement = "";
		replacement += kmer[posInKmer];
		if (errorType == ErrorType::DEL_OF_A) {
			replacement += "A";
		} else if (errorType == ErrorType::DEL_OF_C) {
			replacement += "C";
		} else if (errorType == ErrorType::DEL_OF_G) {
			replacement += "G";
		} else if (errorType == ErrorType::DEL_OF_T) {
			replacement += "T";
		} else {
			throw std::runtime_error("Multidel not supported yet!");
		}
	}
	return replacement;
}

std::string kmerAfterError(const std::string& kmer, ErrorType error, int posOfError) {
	std::string res = kmer;
	res.replace(posOfError, 1, findReplacement(res, error, posOfError));
	return res;
}

std::vector<std::pair<size_t, uint8_t> > rankPotentialErrorPositions(const std::vector<uint8_t>& badKmerCov) {
	std::vector<std::pair<size_t, uint8_t> > res(badKmerCov.size());

	for (size_t i = 0; i < badKmerCov.size(); ++i) {
		res[i] = std::make_pair(i, badKmerCov[i]);
	}
	auto cmp =
			[](const std::pair<size_t, uint8_t>& left, const std::pair<size_t, uint8_t>& right) {return left.second < right.second;};
	std::make_heap(res.begin(), res.end(), cmp);
	/*std::pop_heap(res.begin(), res.end());
	 auto largest = res.back();
	 res.pop_back();*/
	return res;
}

/**
 * For each position, count the number of UNTRUSTED k-mers covered by the position.
 */
std::vector<uint8_t> badKmerCoverage(const std::string& read, CorrectionParameters& params) {
	std::vector<uint8_t> res(read.size());

	size_t i = 0;
	while (i + params.minK < read.size()) {
		size_t k = params.minK;
		std::string kmer = read.substr(i, k);
		KmerType type = classifyKmer(kmer, params.kmerCounter, params.pusm, params.biasUnit,
				params.pathToOriginalReads);

		while (type == KmerType::REPEAT) {
			k += 2;
			if (i + k >= read.size()) {
				break;
			}
			kmer = read.substr(i, k);
			type = classifyKmer(kmer, params.kmerCounter, params.pusm, params.biasUnit, params.pathToOriginalReads);
		}

		if (type == KmerType::UNTRUSTED) {
			for (size_t j = i; j < i + k; ++j) {
				res[j]++;
			}
		}
		i++;
	}

	return res;
}

bool readIsPerfect(const std::string& read, CorrectionParameters& params, size_t from, size_t to) {
	bool perfect = true;
	size_t i = from;

	if (from > read.size() || to > read.size()) {
		throw std::runtime_error(
				"from: " + std::to_string(from) + ", to: " + std::to_string(to) + ", read size: "
						+ std::to_string(read.size()));
	}
	while (i + params.minK <= to) {
		size_t k = params.minK;
		std::string kmer = read.substr(i, k);
		KmerType type = classifyKmer(kmer, params.kmerCounter, params.pusm, params.biasUnit,
				params.pathToOriginalReads);

		while (type == KmerType::REPEAT) {
			k += 2;
			if (i + k >= read.size()) {
				break;
			}
			kmer = read.substr(i, k);
			type = classifyKmer(kmer, params.kmerCounter, params.pusm, params.biasUnit, params.pathToOriginalReads);
		}

		if (type == KmerType::UNTRUSTED) {
			perfect = false;
			break;
		}
		i++;
	}
	return perfect;
}

/*
 * Compute the number of UNTRUSTED k-mers covering a read.
 */
uint8_t numUntrustedKmers(const std::string& read, size_t minK, Matcher& kmerCounter,
		PerfectUniformSequencingModel& pusm, coverage::CoverageBiasUnitMulti& biasUnit,
		const std::string& pathToOriginalReads) {
	uint8_t num = 0;
	size_t i = 0;
	while (i + minK < read.size()) {
		size_t k = minK;
		std::string kmer = read.substr(i, k);
		KmerType type = classifyKmer(kmer, kmerCounter, pusm, biasUnit, pathToOriginalReads);

		while (type == KmerType::REPEAT) {
			k += 2;
			if (i + k >= read.size()) {
				break;
			}
			kmer = read.substr(i, k);
			type = classifyKmer(kmer, kmerCounter, pusm, biasUnit, pathToOriginalReads);
		}

		if (type == KmerType::UNTRUSTED) {
			num++;
		}
		i++;
	}
	return num;
}

/*
 * Compute the number of UNTRUSTED k-mers covering a read[start..end] sequence.
 */
uint8_t numUntrustedKmers(const std::string& read, size_t start, size_t end, CorrectionParameters& params) {
	uint8_t num = 0;
	size_t i = start;
	while (i + params.minK <= end) {
		size_t k = params.minK;
		std::string kmer = read.substr(i, k);
		KmerType type = classifyKmer(kmer, params.kmerCounter, params.pusm, params.biasUnit,
				params.pathToOriginalReads);

		while (type == KmerType::REPEAT) {
			k += 2;
			if (i + k >= read.size()) {
				break;
			}
			kmer = read.substr(i, k);
			type = classifyKmer(kmer, params.kmerCounter, params.pusm, params.biasUnit, params.pathToOriginalReads);
		}

		if (type == KmerType::UNTRUSTED) {
			num++;
		}
		i++;
	}
	return num;
}

size_t findSmallestNonrepetitive(const std::string& str, size_t pos, CorrectionParameters& params) {
	KmerType type = KmerType::REPEAT;
	for (size_t i = params.minK; i + pos < str.size(); i += 2) {
		type = params.classifier.classifyKmer(str.substr(pos, i));
		if (type != KmerType::REPEAT) {
			return i;
		}
	}
	return std::numeric_limits<size_t>::max();
}

Read correctRead_none(const io::Read& read, CorrectionParameters& params) {
	(void) params;
	return read;
}

void performSimpleCorrections(const io::Read& read, CorrectionParameters& params) {
	io::Read correctedRead(read);
	for (size_t i = 0; i < correctedRead.seq.size() - params.minK; ++i) {
		size_t k = findSmallestNonrepetitive(correctedRead.seq, i, params);
		if (k == std::numeric_limits<size_t>::max()) {
			break;
		}
		KmerType type = params.classifier.classifyKmer(correctedRead.seq.substr(i, k));
	}

	throw std::runtime_error("not implemented yet");
}

/*
 * For all error types (including no error), check how many k-mers would still be untrusted... then take the one with the lowest number of untrusted k-mers remaining
 */
bool tryFixingPosition(io::Read& read, size_t pos, CorrectionParameters& params) {
	std::pair<size_t, size_t> affectedArea = affectedReadArea(pos, read.seq.size(), params.minK);
	uint8_t untrustedNormal = numUntrustedKmers(read.seq, affectedArea.first, affectedArea.second, params);
	uint8_t lowestCount = untrustedNormal;
	ErrorType bestType = ErrorType::CORRECT;

	for (ErrorType type : ErrorOnlyTypeIterator()) {
		if (params.correctSingleIndels == false && (isGapErrorType(type) || type == ErrorType::INSERTION)) {
			continue;
		}
		if ((pos == 0 || pos == read.seq.size() - 1) && (isGapErrorType(type) || type == ErrorType::INSERTION)) {
			continue;
		}
		if (params.correctMultidels == false && type == ErrorType::MULTIDEL) {
			continue;
		}

		if (type == ErrorType::SUB_OF_A && read.seq[pos] == 'A') {
			continue;
		}
		if (type == ErrorType::SUB_OF_C && read.seq[pos] == 'C') {
			continue;
		}
		if (type == ErrorType::SUB_OF_G && read.seq[pos] == 'G') {
			continue;
		}
		if (type == ErrorType::SUB_OF_T && read.seq[pos] == 'T') {
			continue;
		}

		std::string readAfterError = kmerAfterError(read.seq, type, pos);
		if (readIsPerfect(readAfterError, params, affectedArea.first, affectedArea.second) && lowestCount > 0) {
			lowestCount = 0;
			bestType = type;
			break;
		}

		/*uint8_t untrustedCount = numUntrustedKmers(readAfterError, affectedArea.first, affectedArea.second, params);
		 if (untrustedCount < lowestCount) {
		 lowestCount = untrustedCount;
		 bestType = type;
		 }*/
	}

	if (bestType != ErrorType::CORRECT /*&& lowestCount == 0*/) {
		read.seq = kmerAfterError(read.seq, bestType, pos);

		//std::cout << errorTypeToString(bestType) << " at " << pos << "\n";
	}

	return (bestType != ErrorType::CORRECT);
}

/*
 * For all error types (including no error), check how many k-mers would still be untrusted... then take the one with the lowest number of untrusted k-mers remaining
 */
bool tryFixingKmer(io::Read& read, size_t pos, size_t k, CorrectionParameters& params) {
	ErrorType bestType = ErrorType::CORRECT;

	std::string kmer = read.seq.substr(pos, k);
	size_t errorPos = 0;

	for (size_t i = 0; i < k; ++i) {
		if (bestType != ErrorType::CORRECT) {
			break;
		}

		for (ErrorType type : ErrorOnlyTypeIterator()) {
			if (params.correctSingleIndels == false && (isGapErrorType(type) || type == ErrorType::INSERTION)) {
				continue;
			}
			if ((i == 0 || i == kmer.size() - 1) && (isGapErrorType(type) || type == ErrorType::INSERTION)) {
				continue;
			}
			if (params.correctMultidels == false && type == ErrorType::MULTIDEL) {
				continue;
			}

			if (type == ErrorType::SUB_OF_A && kmer[i] == 'A') {
				continue;
			}
			if (type == ErrorType::SUB_OF_C && kmer[i] == 'C') {
				continue;
			}
			if (type == ErrorType::SUB_OF_G && kmer[i] == 'G') {
				continue;
			}
			if (type == ErrorType::SUB_OF_T && kmer[i] == 'T') {
				continue;
			}

			std::string kmerNew = kmerAfterError(kmer, type, i);
			if (params.classifier.classifyKmer(kmerNew) != KmerType::UNTRUSTED) {
				bestType = type;
				errorPos = i;
				break;
			}
		}
	}

	if (bestType != ErrorType::CORRECT) {
		read.seq = kmerAfterError(read.seq, bestType, pos + errorPos);
		//std::cout << errorTypeToString(bestType) << " at " << pos << "\n";
	}

	return (bestType != ErrorType::CORRECT);
}

Read correctRead_simple_kmer(const io::Read& read, CorrectionParameters& params) {
	io::Read correctedRead(read);
	bool changed = true;
	while (changed) {
		std::vector<uint8_t> cov = badKmerCoverage(correctedRead.seq, params);
		std::vector<std::pair<size_t, uint8_t> > ranking = rankPotentialErrorPositions(cov);
		changed = false;

		while (ranking.size() > 0) {
			std::pair<size_t, uint8_t> best = ranking.front();
			ranking.pop_back();

			if (best.second > 0) {
				changed = tryFixingPosition(correctedRead, best.first, params);
				if (changed) {
					cov = badKmerCoverage(correctedRead.seq, params);
					ranking = rankPotentialErrorPositions(cov);
					break;
				}
			}
		}

	}
	return correctedRead;
}

Read correctRead_adaptive_kmer(const io::Read& read, CorrectionParameters& params) {
	/*
	 * TODO:
	 * - use sliding window of fixed size? Or just... increase k-mer size as long as k-mer is repetitive?
	 * - check whether highly similar k-mers are also trusted/ other errors are also highly supported
	 * --> if so, further increase k-mer size
	 * - find a good way to deal with conflicting correction candidates from different iterations (this was missing in thesis)
	 * - find a better minimum size for k than 1
	 */
	io::Read correctedRead(read);
	size_t pos = 0;
	while (pos + params.minK < correctedRead.seq.size()) { // loop over starting position
		//find smallest nonrepetitive k-mer size
		size_t k = findSmallestNonrepetitive(correctedRead.seq, pos, params);
		if (k == std::numeric_limits<size_t>::max()) {
			break;
		}
		KmerType type = params.classifier.classifyKmer(correctedRead.seq.substr(pos, k));
		if (type == KmerType::UNTRUSTED) {
			tryFixingKmer(correctedRead, pos, k, params);
		}
		pos += 2;
	}

	return correctedRead;
}

Read correctRead_suffix_tree(const io::Read& read, CorrectionParameters& params) {
	throw std::runtime_error("not implemented yet");
}

Read correctRead_full_msa(const io::Read& read, CorrectionParameters& params) {
	throw std::runtime_error("not implemented yet");
}

Read correctRead_partial_msa(const io::Read& read, CorrectionParameters& params) {
	throw std::runtime_error("not implemented yet");
}

Read correctRead(const Read& read, CorrectionAlgorithm algo, CorrectionParameters& params) {
	switch (algo) {
	case CorrectionAlgorithm::NONE:
		return correctRead_none(read, params);
		break;
	case CorrectionAlgorithm::SIMPLE_KMER:
		return correctRead_simple_kmer(read, params);
		break;
	case CorrectionAlgorithm::ADAPTIVE_KMER:
		return correctRead_adaptive_kmer(read, params);
		break;
	case CorrectionAlgorithm::FULL_MSA:
		return correctRead_full_msa(read, params);
		break;
	case CorrectionAlgorithm::PARTIAL_MSA:
		return correctRead_partial_msa(read, params);
		break;
	case CorrectionAlgorithm::SUFFIX_TREE:
		return correctRead_suffix_tree(read, params);
		break;
	default:
		throw std::runtime_error("Unknown correction algorithm");
	}
}

void correctReads(const std::string& pathToOriginalReads, CorrectionAlgorithm algo, Matcher& kmerCounter,
		PerfectUniformSequencingModel& pusm, coverage::CoverageBiasUnitMulti& biasUnit, const std::string& outputPath,
		size_t kmerSize) {
	//std::cout << "count of " << "ACTCACTGTTT:" << kmerCounter.countKmer("ACTCACTGTTT");

	bool correctSingleIndels = true;
	bool correctMultidels = false;

	std::cout << "correctSIngleIndels: " << correctSingleIndels << "\n";

	io::ReadOutput printer;
	printer.createFile(outputPath);

	io::ReadInput reader(pathToOriginalReads);
	HashClassifier classifier(kmerCounter, pusm, biasUnit, pathToOriginalReads);
	biasUnit.preprocess(kmerSize, pathToOriginalReads, kmerCounter, pusm);
	biasUnit.plotMedianCoverageBiases(outputPath + "_biases");
	biasUnit.preprocess(kmerSize+2, pathToOriginalReads, kmerCounter, pusm);
		biasUnit.plotMedianCoverageBiases(outputPath + "_biases");
		biasUnit.preprocess(kmerSize+4, pathToOriginalReads, kmerCounter, pusm);
			biasUnit.plotMedianCoverageBiases(outputPath + "_biases");
			biasUnit.preprocess(kmerSize+6, pathToOriginalReads, kmerCounter, pusm);
						biasUnit.plotMedianCoverageBiases(outputPath + "_biases");
	CorrectionParameters params(kmerSize, kmerCounter, pusm, biasUnit, pathToOriginalReads, classifier,
			correctSingleIndels, correctMultidels);

#pragma omp parallel
	{
#pragma omp single
		{
			while (reader.hasNext()) {
				io::Read uncorrected = reader.readNext(true, false, true);
#pragma omp task shared(printer, params) firstprivate(uncorrected)
				{
					io::Read corrected = correctRead(uncorrected, algo, params);
#pragma omp critical
					printer.write(corrected);
				}
			}
		}
#pragma omp taskwait
	}
}

} // end of namespace seq_correct::correction
} // end of namespace seq_correct
