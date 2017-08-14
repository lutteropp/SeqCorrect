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

#include <unordered_map>
#include <utility>
#include "../util/enums.hpp"

namespace seq_correct {
namespace eval {

using namespace util;

// row in confusion matrix : true error type
// column in confusion matrix: predicted error type
class ErrorEvaluationData {
public:
	ErrorEvaluationData();
	size_t truePositives(ErrorType type) const;
	size_t falsePositives(ErrorType type) const;
	size_t trueNegatives(ErrorType type) const;
	size_t falseNegatives(ErrorType type) const;
	size_t getEntry(ErrorType trueType, ErrorType predictedType) const;
	size_t sumAllBaseEntries() const;
	size_t sumAllGapEntries() const;
	size_t sumTruth(ErrorType type) const; // row sum
	size_t sumPredicted(ErrorType type) const; // column sum
	void update(ErrorType trueType, ErrorType predictedType);
private:
	std::unordered_map<std::pair<ErrorType, ErrorType>, size_t, EnumClassPairHash> baseConfusionMatrix;
	std::unordered_map<std::pair<ErrorType, ErrorType>, size_t, EnumClassPairHash> gapConfusionMatrix;
};

class KmerEvaluationData {
public:
	KmerEvaluationData();
	size_t truePositives(KmerType type) const;
	size_t falsePositives(KmerType type) const;
	size_t trueNegatives(KmerType type) const;
	size_t falseNegatives(KmerType type) const;
	size_t getEntry(KmerType trueType, KmerType predictedType) const;
	size_t sumAllEntries() const;
	size_t sumTruth(KmerType type) const; // row sum
	size_t sumPredicted(KmerType type) const; // column sum
	void update(KmerType trueType, KmerType predictedType);
private:
	std::unordered_map<std::pair<KmerType, KmerType>, size_t, EnumClassPairHash> kmerConfusionMatrix;
};

} // end of namespace seq_correct::eval
} // end of namespace seq_correct
