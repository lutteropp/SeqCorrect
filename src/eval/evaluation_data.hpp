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
#include "../util/error_type.hpp"

namespace seq_correct {
namespace eval {

using namespace util;

// row in confusion matrix : true error type
// column in confusion matrix: predicted error type
struct EvaluationData {
	EvaluationData();
	size_t truePositives(ErrorType type) const;
	size_t falsePositives(ErrorType type) const;
	size_t trueNegatives(ErrorType type) const;
	size_t falseNegatives(ErrorType type) const;

	std::unordered_map<std::pair<ErrorType, ErrorType>, size_t, EnumClassPairHash> baseConfusionMatrix;
	std::unordered_map<std::pair<ErrorType, ErrorType>, size_t, EnumClassPairHash> gapConfusionMatrix;
};

} // end of namespace seq_correct::eval
} // end of namespace seq_correct
