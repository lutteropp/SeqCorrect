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
#include <unordered_map>
#include <utility>
#include "../util/error_type.hpp"

namespace seq_correct {
namespace eval {

using namespace util;

struct EvaluationData {
	EvaluationData() {
		for (ErrorType e1 : AllBaseTypeIterator()) {
			for (ErrorType e2 : AllBaseTypeIterator()) {
				baseConfusionMatrix[std::make_pair(e1, e2)] = 0;
			}
		}
		for (ErrorType e1 : AllGapTypeIterator()) {
			for (ErrorType e2 : AllGapTypeIterator()) {
				gapConfusionMatrix[std::make_pair(e1, e2)] = 0;
			}
		}
		for (ErrorType e : AllErrorTypeIterator() ) {
			truePositives[e] = 0;
			falsePositives[e] = 0;
			trueNegatives[e] = 0;
			truePositives[e] = 0;
		}
	}
	std::unordered_map<ErrorType, size_t> truePositives;
	std::unordered_map<ErrorType, size_t> falsePositives;
	std::unordered_map<ErrorType, size_t> trueNegatives;
	std::unordered_map<ErrorType, size_t> falseNegatives;
	std::unordered_map<std::pair<ErrorType, ErrorType>, size_t > baseConfusionMatrix;
	std::unordered_map<std::pair<ErrorType, ErrorType>, size_t > gapConfusionMatrix;
};

double computeAccuracy(ErrorType type, const EvaluationData& data);
double computePrecision(ErrorType type, const EvaluationData& data);
double computeRecall(ErrorType type, const EvaluationData& data);
double computeSpecificity(ErrorType type, const EvaluationData& data);
double computeF1Score(ErrorType type, const EvaluationData& data);
double computeBaseNMIScore(const EvaluationData& data);
double computeGapNMIScore(const EvaluationData& data);
double computeUnweightedAverageBaseF1Score(const EvaluationData& data);
double computeUnweightedAverageGapF1Score(const EvaluationData& data);

EvaluationData evaluateCorrections(const std::string& originalReadsFilepath, const std::string& correctedReadsFilepath,
		const std::string& genomeFilepath);

EvaluationData evaluateCorrections(const std::string& alignmentFilepath, const std::string& correctedReadsFilepath);

} // end of namespace seq_correct::eval
} // end of namespace seq_correct
