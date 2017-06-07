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

#include <utility>

namespace seq_correct {
namespace util {

// code from http://stackoverflow.com/questions/261963/how-can-i-iterate-over-an-enum
template<typename C, C beginVal, C endVal>
class Iterator {
	typedef typename std::underlying_type<C>::type val_t;
	int val;
public:
	Iterator(const C & f) :
			val(static_cast<val_t>(f)) {
	}
	Iterator() :
			val(static_cast<val_t>(beginVal)) {
	}
	Iterator operator++() {
		++val;
		return *this;
	}
	C operator*() {
		return static_cast<C>(val);
	}
	Iterator begin() {
		return *this;
	} //default ctor is good
	Iterator end() {
		static const Iterator endIter = ++Iterator(endVal); // cache it
		return endIter;
	}
	bool operator!=(const Iterator& i) {
		return val != i.val;
	}
};

//modified from https://stackoverflow.com/questions/18837857/cant-use-enum-class-as-unordered-map-key
//  and https://stackoverflow.com/questions/11376163/custom-hash-function-for-pair-of-enum-values-used-as-unordered-map-key
struct EnumClassPairHash
{
    template <typename T>
    std::size_t operator()(std::pair<T, T> t) const
    {
        return static_cast<std::size_t>(t.first) ^ static_cast<std::size_t>(t.second);
    }
};

enum class ErrorType {
	CORRECT, INSERTION, SUB_OF_A, SUB_OF_C, SUB_OF_G, SUB_OF_T, DEL_OF_A, DEL_OF_C, DEL_OF_G, DEL_OF_T, MULTIDEL, NODEL
};
typedef Iterator<ErrorType, ErrorType::CORRECT, ErrorType::NODEL> AllErrorTypeIterator;
typedef Iterator<ErrorType, ErrorType::INSERTION, ErrorType::MULTIDEL> ErrorOnlyTypeIterator;
typedef Iterator<ErrorType, ErrorType::CORRECT, ErrorType::SUB_OF_T> BaseTypeIterator;
typedef Iterator<ErrorType, ErrorType::DEL_OF_A, ErrorType::NODEL> GapTypeIterator;

inline bool isBaseErrorType(ErrorType type) {
	for (ErrorType e : BaseTypeIterator()) {
		if (e == type) {
			return true;
		}
	}
	return false;
}

inline bool isGapErrorType(ErrorType type) {
	for (ErrorType e : GapTypeIterator()) {
		if (e == type) {
			return true;
		}
	}
	return false;
}

} // end of namespace seq_correct::util
} // end of namespace seq_correct
