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

#include "util.hpp"

#include <stdexcept>
#include <algorithm>
#include <iostream>

namespace seq_correct {
namespace util {

io::Read reverseComplementRead(const io::Read& read) {
	io::Read rcRead;
	rcRead.name = read.name;
	rcRead.seq = reverseComplementString(read.seq);
	rcRead.qual = read.qual;
	std::reverse(rcRead.qual.begin(), rcRead.qual.end());
	return rcRead;
}

double gcContent(const std::string& kmer) {
	size_t gcCount = std::count_if(kmer.begin(), kmer.end(), [](char c) {return (c=='C') || (c=='G');});
	return gcCount / (double) kmer.size();
}

std::string kmerAfterError(const std::string& kmer, size_t pos, ErrorType type) {
	if (pos >= kmer.size()) {
		throw std::runtime_error("pos >= kmer.size()");
	}
	std::string newKmer;
	switch (type) {
	case ErrorType::CORRECT:
		newKmer = kmer;
		break;
	case ErrorType::INSERTION:
		newKmer = kmer.substr(0, pos) + kmer.substr(pos + 1, std::string::npos);
		break;
	case ErrorType::SUB_OF_A:
		newKmer = kmer;
		newKmer[pos] = 'A';
		break;
	case ErrorType::SUB_OF_C:
		newKmer = kmer;
		newKmer[pos] = 'C';
		break;
	case ErrorType::SUB_OF_G:
		newKmer = kmer;
		newKmer[pos] = 'G';
		break;
	case ErrorType::SUB_OF_T:
		newKmer = kmer;
		newKmer[pos] = 'T';
		break;
	case ErrorType::DEL_OF_A:
		newKmer = kmer.substr(0, pos) + "A" + kmer.substr(pos, std::string::npos);
		break;
	case ErrorType::DEL_OF_C:
		newKmer = kmer.substr(0, pos) + "C" + kmer.substr(pos, std::string::npos);
		break;
	case ErrorType::DEL_OF_G:
		newKmer = kmer.substr(0, pos) + "G" + kmer.substr(pos, std::string::npos);
		break;
	case ErrorType::DEL_OF_T:
		newKmer = kmer.substr(0, pos) + "T" + kmer.substr(pos, std::string::npos);
		break;
	case ErrorType::MULTIDEL:
		if (kmer[pos] == '_') {
			throw std::runtime_error("multidel after a multidel");
		}
		newKmer = kmer.substr(0, pos) + "_" + kmer.substr(pos, std::string::npos);
		break;
	case ErrorType::NODEL:
		newKmer = kmer;
		break;
	}
	return newKmer;
}

} // end of namespace seq_correct::util
} // end of namespace seq_correct
