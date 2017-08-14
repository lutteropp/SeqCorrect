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

#include <cmath>
#include "metrics.hpp"

namespace seq_correct {
namespace eval {

size_t truePositives(ErrorType type, const ErrorEvaluationData& data) {
	return data.truePositives(type);
}

size_t trueNegatives(ErrorType type, const ErrorEvaluationData& data) {
	return data.trueNegatives(type);
}

size_t falsePositives(ErrorType type, const ErrorEvaluationData& data) {
	return data.falsePositives(type);
}

size_t falseNegatives(ErrorType type, const ErrorEvaluationData& data) {
	return data.falseNegatives(type);
}

size_t truePositives(KmerType type, const KmerEvaluationData& data) {
	return data.truePositives(type);
}

size_t trueNegatives(KmerType type, const KmerEvaluationData& data) {
	return data.trueNegatives(type);
}

size_t falsePositives(KmerType type, const KmerEvaluationData& data) {
	return data.falsePositives(type);
}

size_t falseNegatives(KmerType type, const KmerEvaluationData& data) {
	return data.falseNegatives(type);
}

double computeAccuracy(ErrorType type, const ErrorEvaluationData& data) {
	return (data.truePositives(type) + data.falsePositives(type))
			/ (double) (data.truePositives(type) + data.falsePositives(type) + data.falsePositives(type)
					+ data.falseNegatives(type));
}

double computeAccuracy(KmerType type, const KmerEvaluationData& data) {
	return (data.truePositives(type) + data.falsePositives(type))
			/ (double) (data.truePositives(type) + data.falsePositives(type) + data.falsePositives(type)
					+ data.falseNegatives(type));
}

double computePrecision(ErrorType type, const ErrorEvaluationData& data) {
	return data.truePositives(type) / (double) (data.truePositives(type) + data.falsePositives(type));
}

double computePrecision(KmerType type, const KmerEvaluationData& data) {
	return data.truePositives(type) / (double) (data.truePositives(type) + data.falsePositives(type));
}


double computeRecall(ErrorType type, const ErrorEvaluationData& data) {
	return data.truePositives(type) / (double) (data.truePositives(type) + data.falseNegatives(type));
}

double computeRecall(KmerType type, const KmerEvaluationData& data) {
	return data.truePositives(type) / (double) (data.truePositives(type) + data.falseNegatives(type));
}

double computeSensitivity(ErrorType type, const ErrorEvaluationData& data) {
	return computeRecall(type, data);
}

double computeSensitivity(KmerType type, const KmerEvaluationData& data) {
	return computeRecall(type, data);
}

double computeGain(ErrorType type, const ErrorEvaluationData& data) {
	size_t tp = data.truePositives(type);
	size_t fp = data.falsePositives(type);
	size_t fn = data.falseNegatives(type);
	return (tp - fp) / (double) (tp + fn);
}

double computeGain(KmerType type, const KmerEvaluationData& data) {
	size_t tp = data.truePositives(type);
	size_t fp = data.falsePositives(type);
	size_t fn = data.falseNegatives(type);
	return (tp - fp) / (double) (tp + fn);
}

double computeSpecificity(ErrorType type, const ErrorEvaluationData& data) {
	return data.trueNegatives(type) / (double) (data.trueNegatives(type) + data.falseNegatives(type));
}

double computeSpecificity(KmerType type, const KmerEvaluationData& data) {
	return data.trueNegatives(type) / (double) (data.trueNegatives(type) + data.falseNegatives(type));
}

double computeF1Score(ErrorType type, const ErrorEvaluationData& data) {
	double precision = computePrecision(type, data);
	double recall = computeRecall(type, data);
	return 2 * (precision * recall) / (precision + recall);
}

double computeF1Score(KmerType type, const KmerEvaluationData& data) {
	double precision = computePrecision(type, data);
	double recall = computeRecall(type, data);
	return 2 * (precision * recall) / (precision + recall);
}


size_t trueTotal(ErrorType type, const ErrorEvaluationData& data) {
	return data.truePositives(type) + data.falseNegatives(type);
}

size_t trueTotal(KmerType type, const KmerEvaluationData& data) {
	return data.truePositives(type) + data.falseNegatives(type);
}

double computeBaseEntropyTruth(const ErrorEvaluationData& data) {
	double sum = 0;
	size_t N = data.sumAllBaseEntries();
	for (ErrorType e : BaseTypeIterator()) {
		double val = data.sumTruth(e) / (double) N;
		if (val != 0) {
			sum += val * std::log2(val);
		}
	}
	sum *= -1;
	return sum;
}

double computeBaseEntropyPredicted(const ErrorEvaluationData& data) {
	double sum = 0;
	size_t N = data.sumAllBaseEntries();
	for (ErrorType e : BaseTypeIterator()) {
		double val = data.sumPredicted(e) / (double) N;
		if (val != 0) {
			sum += val * std::log2(val);
		}
	}
	sum *= -1;
	return sum;
}

double computeGapEntropyTruth(const ErrorEvaluationData& data) {
	double sum = 0;
	size_t N = data.sumAllGapEntries();
	for (ErrorType e : GapTypeIterator()) {
		double val = data.sumTruth(e) / (double) N;
		if (val != 0) {
			sum += val * std::log2(val);
		}
	}
	sum *= -1;
	return sum;
}

double computeGapEntropyPredicted(const ErrorEvaluationData& data) {
	double sum = 0;
	size_t N = data.sumAllGapEntries();
	for (ErrorType e : GapTypeIterator()) {
		double val = data.sumPredicted(e) / (double) N;
		if (val != 0) {
			sum += val * std::log2(val);
		}
	}
	sum *= -1;
	return sum;
}

double computeEntropyTruth(const KmerEvaluationData& data) {
	double sum = 0;
	size_t N = data.sumAllEntries();
	for (KmerType e : KmerTypeIterator()) {
		double val = data.sumTruth(e) / (double) N;
		if (val != 0) {
			sum += val * std::log2(val);
		}
	}
	sum *= -1;
	return sum;
}

double computeEntropyPredicted(const KmerEvaluationData& data) {
	double sum = 0;
	size_t N = data.sumAllEntries();
	for (KmerType e : KmerTypeIterator()) {
		double val = data.sumPredicted(e) / (double) N;
		if (val != 0) {
			sum += val * std::log2(val);
		}
	}
	sum *= -1;
	return sum;
}

double computeGapMutualInformation(const ErrorEvaluationData& data) {
	double mi = 0;
	size_t N = data.sumAllGapEntries();
	size_t NSquare = N * N;
	for (ErrorType e : GapTypeIterator()) {
		for (ErrorType f : GapTypeIterator()) {
			double part1 = data.getEntry(e, f) / (double) N;
			double part2 = data.sumTruth(e) * data.sumPredicted(f) / (double) (NSquare);
			if (part1 != 0 && part2 != 0) {
				mi += part1 * std::log2(part1 / part2);
			}
		}
	}
	return mi;
}

double computeBaseMutualInformation(const ErrorEvaluationData& data) {
	double mi = 0;
	size_t N = data.sumAllBaseEntries();
	size_t NSquare = N * N;
	for (ErrorType e : BaseTypeIterator()) {
		for (ErrorType f : BaseTypeIterator()) {
			double part1 = data.getEntry(e, f) / (double) N;
			double part2 = data.sumTruth(e) * data.sumPredicted(f) / (double) (NSquare);
			if (part1 != 0 && part2 != 0) {
				mi += part1 * std::log2(part1 / part2);
			}
		}
	}
	return mi;
}

double computeMutualInformation(const KmerEvaluationData& data) {
	double mi = 0;
	size_t N = data.sumAllEntries();
	size_t NSquare = N * N;
	for (KmerType e : KmerTypeIterator()) {
		for (KmerType f : KmerTypeIterator()) {
			double part1 = data.getEntry(e, f) / (double) N;
			double part2 = data.sumTruth(e) * data.sumPredicted(f) / (double) (NSquare);
			if (part1 != 0 && part2 != 0) {
				mi += part1 * std::log2(part1 / part2);
			}
		}
	}
	return mi;
}

// Using the NMI_max formula I(U,V) / max(H(U), H(V))
double computeBaseNMIScore(const ErrorEvaluationData& data) {
	return computeBaseMutualInformation(data)
			/ std::max(computeBaseEntropyTruth(data), computeBaseEntropyPredicted(data));
}

// Using the NMI_max formula I(U,V) / max(H(U), H(V))
double computeGapNMIScore(const ErrorEvaluationData& data) {
	return computeGapMutualInformation(data) / std::max(computeGapEntropyTruth(data), computeGapEntropyPredicted(data));
}

// Using the NMI_max formula I(U,V) / max(H(U), H(V))
double computeNMIScore(const KmerEvaluationData& data) {
	return computeMutualInformation(data)
			/ std::max(computeEntropyTruth(data), computeEntropyPredicted(data));
}

double computeUnweightedAverageBaseF1Score(const ErrorEvaluationData& data) {
	double fscoreSum = 0.0;
	int numValid = 0;
	double fscore = computeF1Score(ErrorType::INSERTION, data);
	if (fscore == fscore) { // check for NaN
		fscoreSum += fscore;
		numValid++;
	}
	fscore = computeF1Score(ErrorType::SUB_OF_A, data);
	if (fscore == fscore) { // check for NaN
		fscoreSum += fscore;
		numValid++;
	}
	fscore = computeF1Score(ErrorType::SUB_OF_C, data);
	if (fscore == fscore) { // check for NaN
		fscoreSum += fscore;
		numValid++;
	}
	fscore = computeF1Score(ErrorType::SUB_OF_G, data);
	if (fscore == fscore) { // check for NaN
		fscoreSum += fscore;
		numValid++;
	}
	fscore = computeF1Score(ErrorType::SUB_OF_T, data);
	if (fscore == fscore) { // check for NaN
		fscoreSum += fscore;
		numValid++;
	}
	return fscoreSum / (double) numValid;
}

double computeUnweightedAverageGapF1Score(const ErrorEvaluationData& data) {
	double fscoreSum = 0.0;
	int numValid = 0;
	double fscore = computeF1Score(ErrorType::MULTIDEL, data);
	if (fscore == fscore) { // check for NaN
		fscoreSum += fscore;
		numValid++;
	}
	fscore = computeF1Score(ErrorType::DEL_OF_A, data);
	if (fscore == fscore) { // check for NaN
		fscoreSum += fscore;
		numValid++;
	}
	fscore = computeF1Score(ErrorType::DEL_OF_C, data);
	if (fscore == fscore) { // check for NaN
		fscoreSum += fscore;
		numValid++;
	}
	fscore = computeF1Score(ErrorType::DEL_OF_G, data);
	if (fscore == fscore) { // check for NaN
		fscoreSum += fscore;
		numValid++;
	}
	fscore = computeF1Score(ErrorType::DEL_OF_T, data);
	if (fscore == fscore) { // check for NaN
		fscoreSum += fscore;
		numValid++;
	}
	return fscoreSum / (double) numValid;
}

double computeUnweightedAverageF1Score(const KmerEvaluationData& data) {
	double fscoreSum = 0.0;
	int numValid = 0;
	double fscore = computeF1Score(KmerType::UNTRUSTED, data);
	if (fscore == fscore) { // check for NaN
		fscoreSum += fscore;
		numValid++;
	}
	fscore = computeF1Score(KmerType::UNIQUE, data);
	if (fscore == fscore) { // check for NaN
		fscoreSum += fscore;
		numValid++;
	}
	fscore = computeF1Score(KmerType::REPEAT, data);
	if (fscore == fscore) { // check for NaN
		fscoreSum += fscore;
		numValid++;
	}
	return fscoreSum / (double) numValid;
}

} // end of namespace seq_correct::eval
} // end of namespace seq_correct
