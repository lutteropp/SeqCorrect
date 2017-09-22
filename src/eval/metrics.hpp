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

#include "../util/enums.hpp"
#include "evaluation_data.hpp"

namespace seq_correct {
namespace eval {

using namespace util;

size_t truePositives(ErrorType type, const ErrorEvaluationData& data);
size_t trueNegatives(ErrorType type, const ErrorEvaluationData& data);
size_t falsePositives(ErrorType type, const ErrorEvaluationData& data);
size_t falseNegatives(ErrorType type, const ErrorEvaluationData& data);

double computeAccuracy(ErrorType type, const ErrorEvaluationData& data);
double computePrecision(ErrorType type, const ErrorEvaluationData& data);
double computeRecall(ErrorType type, const ErrorEvaluationData& data);
double computeSensitivity(ErrorType type, const ErrorEvaluationData& data);
double computeGain(ErrorType type, const ErrorEvaluationData& data);
double computeSpecificity(ErrorType type, const ErrorEvaluationData& data);
double computeF1Score(ErrorType type, const ErrorEvaluationData& data);
double computeUnweightedAverageBaseF1Score(const ErrorEvaluationData& data);
double computeUnweightedAverageGapF1Score(const ErrorEvaluationData& data);
double computeBaseNMIScore(const ErrorEvaluationData& data);
double computeGapNMIScore(const ErrorEvaluationData& data);

size_t number_confused_errors(const ErrorEvaluationData& data);
size_t number_introduced_errors(const ErrorEvaluationData& data);
size_t number_claimed_errors(const ErrorEvaluationData& data);
size_t number_undiscovered_errors(const ErrorEvaluationData& data);
size_t number_discovered_errors(const ErrorEvaluationData& data);

size_t truePositives(KmerType type, const KmerEvaluationData& data);
size_t trueNegatives(KmerType type, const KmerEvaluationData& data);
size_t falsePositives(KmerType type, const KmerEvaluationData& data);
size_t falseNegatives(KmerType type, const KmerEvaluationData& data);

double computeAccuracy(KmerType type, const KmerEvaluationData& data);
double computePrecision(KmerType type, const KmerEvaluationData& data);
double computeRecall(KmerType type, const KmerEvaluationData& data);
double computeSensitivity(KmerType type, const KmerEvaluationData& data);
double computeGain(KmerType type, const KmerEvaluationData& data);
double computeSpecificity(KmerType type, const KmerEvaluationData& data);
double computeF1Score(KmerType type, const KmerEvaluationData& data);
double computeUnweightedAverageF1Score(const KmerEvaluationData& data);
double computeNMIScore(const KmerEvaluationData& data);

} // end of namespace seq_correct::eval
} // end of namespace seq_correct
