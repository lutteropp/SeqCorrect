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

#include "../external/const_string_ptr.hpp"
#include "../io/sequence_io.hpp"

namespace seq_correct {
namespace util {

	enum class ErrorType {
		INSERTION, SUB_OF_A, SUB_OF_C, SUB_OF_G, SUB_OF_T, DEL_OF_A, DEL_OF_C, DEL_OF_G, DEL_OF_T, MULTIDEL
	};

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

	typedef Iterator<ErrorType, ErrorType::INSERTION, ErrorType::MULTIDEL> ErrorTypeIterator;

	std::string reverseComplementString(const std::string& text);
	io::Read reverseComplementRead(const io::Read& read);
	double gcContent(const std::string& kmer);

	template <typename T>
	size_t countGC(const T& kmer) {
		size_t gcCount = 0;
		for (size_t i = 0; i < kmer.size(); ++i) {
			if (kmer[i] == 'G' || kmer[i] == 'C') {
				gcCount++;
			}
		}
		return gcCount;
	}

	std::string kmerAfterError(const std::string& kmer, size_t pos, ErrorType type);

	enum class GenomeType {
		CIRCULAR, LINEAR
	};

	class Dataset {
	public:
		Dataset(GenomeType genomeType, size_t genomeSize, const std::string& readFilepath);
	private:
		GenomeType genomeType;
		size_t genomeSize;
		std::unordered_map<size_t, size_t> readLengths;
		const std::string& readFilepath;
	};

	class ReferenceDataset: public Dataset {
	public:
		ReferenceDataset(GenomeType genomeType, size_t genomeSize, const std::string& readFilepath,
				const std::string& referenceGenomePath);
	private:
		std::string referenceGenomePath;
		std::string referenceGenome;
	};

} // end of namespace seq_correct::util
} // end of namespace seq_correct
