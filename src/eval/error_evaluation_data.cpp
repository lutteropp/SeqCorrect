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

#include "error_evaluation_data.hpp"

#include <stdexcept>
#include <iostream>

namespace seq_correct {
namespace eval {

using namespace util;

KmerEvaluationData::KmerEvaluationData() {
	for (KmerType k1: KmerTypeIterator()) {
		for (KmerType k2: KmerTypeIterator()) {
			kmerConfusionMatrix[std::make_pair(k1, k2)] = 0;
		}
	}
}

ErrorEvaluationData::ErrorEvaluationData() {
	for (ErrorType e1 : BaseTypeIterator()) {
		for (ErrorType e2 : BaseTypeIterator()) {
			baseConfusionMatrix[std::make_pair(e1, e2)] = 0;
		}
	}
	for (ErrorType e1 : GapTypeIterator()) {
		for (ErrorType e2 : GapTypeIterator()) {
			gapConfusionMatrix[std::make_pair(e1, e2)] = 0;
		}
	}
}

size_t KmerEvaluationData::truePositives(KmerType type) const {
	return kmerConfusionMatrix.at(std::make_pair(type, type));
}

size_t ErrorEvaluationData::truePositives(ErrorType type) const {
	if (isBaseErrorType(type)) {
		return baseConfusionMatrix.at(std::make_pair(type, type));
	} else {
		return gapConfusionMatrix.at(std::make_pair(type, type));
	}
}

size_t KmerEvaluationData::falsePositives(KmerType type) const {
	size_t count = 0;
	for (KmerType t : KmerTypeIterator()) {
		if (t != type) {
			count += kmerConfusionMatrix.at(std::make_pair(t, type));
		}
	}
	return count;
}

size_t ErrorEvaluationData::falsePositives(ErrorType type) const {
	size_t count = 0;
	if (isBaseErrorType(type)) {
		for (ErrorType e : BaseTypeIterator()) {
			if (e != type) {
				count += baseConfusionMatrix.at(std::make_pair(e, type));
			}
		}
	} else {
		for (ErrorType e : GapTypeIterator()) {
			if (e != type) {
				count += gapConfusionMatrix.at(std::make_pair(e, type));
			}
		}
	}
	return count;
}

size_t KmerEvaluationData::trueNegatives(KmerType type) const {
	size_t count = 0;
	for (KmerType t : KmerTypeIterator()) {
		if (t != type) {
			for (KmerType u : KmerTypeIterator()) {
				if (u != type) {
					count += kmerConfusionMatrix.at(std::make_pair(t, u));
				}
			}
		}
	}
	return count;
}

size_t ErrorEvaluationData::trueNegatives(ErrorType type) const {
	size_t count = 0;
	if (isBaseErrorType(type)) {
		for (ErrorType e : BaseTypeIterator()) {
			if (e != type) {
				for (ErrorType f : BaseTypeIterator()) {
					if (f != type) {
						count += baseConfusionMatrix.at(std::make_pair(e, f));
					}
				}
			}
		}
	} else {
		for (ErrorType e : GapTypeIterator()) {
			if (e != type) {
				for (ErrorType f : GapTypeIterator()) {
					if (f != type) {
						count += gapConfusionMatrix.at(std::make_pair(e, f));
					}
				}
			}
		}
	}
	return count;
}

size_t KmerEvaluationData::falseNegatives(KmerType type) const {
	size_t count = 0;
	for (KmerType t : KmerTypeIterator()) {
		if (t != type) {
			count += kmerConfusionMatrix.at(std::make_pair(type, t));
		}
	}
	return count;
}

size_t ErrorEvaluationData::falseNegatives(ErrorType type) const {
	size_t count = 0;
	if (isBaseErrorType(type)) {
		for (ErrorType e : BaseTypeIterator()) {
			if (e != type) {
				count += baseConfusionMatrix.at(std::make_pair(type, e));
			}
		}
	} else {
		for (ErrorType e : GapTypeIterator()) {
			if (e != type) {
				count += gapConfusionMatrix.at(std::make_pair(type, e));
			}
		}
	}
	return count;
}

size_t KmerEvaluationData::getEntry(KmerType trueType, KmerType predictedType) const {
	return kmerConfusionMatrix.at(std::make_pair(trueType, predictedType));
}

size_t ErrorEvaluationData::getEntry(ErrorType trueType, ErrorType predictedType) const {
	if (isBaseErrorType(trueType) && isBaseErrorType(predictedType)) {
		return baseConfusionMatrix.at(std::make_pair(trueType, predictedType));
	} else if (isGapErrorType(trueType) && isGapErrorType(predictedType)) {
		return gapConfusionMatrix.at(std::make_pair(trueType, predictedType));
	} else {
		std::cout << util::errorTypeToString(trueType) << ":\n";
		if (isBaseErrorType(trueType)) {
			std::cout << "trueType is a base error type\n";
		}
		if (isGapErrorType(trueType)) {
			std::cout << "trueType is a gap error type\n";
		}
		std::cout << util::errorTypeToString(predictedType) << ":\n";
		if (isBaseErrorType(predictedType)) {
			std::cout << "predictedType is a base error type\n";
		}
		if (isGapErrorType(predictedType)) {
			std::cout << "predictedType is a gap error type\n";
		}
		throw std::runtime_error("trueType and predictedType must be both base types or both gap types!");
	}
}

size_t ErrorEvaluationData::sumAllBaseEntries() const {
	size_t count = 0;
	for (ErrorType e : BaseTypeIterator()) {
		for (ErrorType f : BaseTypeIterator()) {
			count += baseConfusionMatrix.at(std::make_pair(e, f));
		}
	}
	return count;
}

size_t ErrorEvaluationData::sumAllGapEntries() const {
	size_t count = 0;
	for (ErrorType e : GapTypeIterator()) {
		for (ErrorType f : GapTypeIterator()) {
			count += gapConfusionMatrix.at(std::make_pair(e, f));
		}
	}
	return count;
}

size_t KmerEvaluationData::sumAllEntries() const {
	size_t count = 0;
	for (KmerType e : KmerTypeIterator()) {
		for (KmerType f : KmerTypeIterator()) {
			count += kmerConfusionMatrix.at(std::make_pair(e, f));
		}
	}
	return count;
}

size_t KmerEvaluationData::sumTruth(KmerType type) const {
	size_t sum = 0;
	for (KmerType t : KmerTypeIterator()) {
		sum += kmerConfusionMatrix.at(std::make_pair(type, t));
	}
	return sum;
}

size_t ErrorEvaluationData::sumTruth(ErrorType type) const {
	size_t sum = 0;
	if (isBaseErrorType(type)) {
		for (ErrorType e : BaseTypeIterator()) {
			sum += baseConfusionMatrix.at(std::make_pair(type, e));
		}
	} else {
		for (ErrorType e : GapTypeIterator()) {
			sum += gapConfusionMatrix.at(std::make_pair(type, e));
		}
	}
	return sum;
}

size_t KmerEvaluationData::sumPredicted(KmerType type) const {
	size_t sum = 0;
	for (KmerType t : KmerTypeIterator()) {
		sum += kmerConfusionMatrix.at(std::make_pair(t, type));
	}
	return sum;
}

size_t ErrorEvaluationData::sumPredicted(ErrorType type) const {
	size_t sum = 0;
	if (isBaseErrorType(type)) {
		for (ErrorType e : BaseTypeIterator()) {
			sum += baseConfusionMatrix.at(std::make_pair(e, type));
		}
	} else {
		for (ErrorType e : GapTypeIterator()) {
			sum += gapConfusionMatrix.at(std::make_pair(e, type));
		}
	}
	return sum;
}

void KmerEvaluationData::update(KmerType trueType, KmerType predictedType) {
#pragma omp critical(updateKmer)
	kmerConfusionMatrix[std::make_pair(trueType, predictedType)]++;
}

void ErrorEvaluationData::update(ErrorType trueType, ErrorType predictedType) {
	if (isBaseErrorType(trueType) && isBaseErrorType(predictedType)) {
#pragma omp critical(updateBase)
		baseConfusionMatrix[std::make_pair(trueType, predictedType)]++;
	} else if (isGapErrorType(trueType) && isGapErrorType(predictedType)) {
#pragma omp critical(updateGap)
		gapConfusionMatrix[std::make_pair(trueType, predictedType)]++;
	} else {
		throw std::runtime_error("trueType and predictedType are of incompatible types");
	}
}

} // end of namespace seq_correct::eval
} // end of namespace seq_correct
