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

namespace seq_correct {
namespace util {

std::string reverseComplementString(const std::string& text) {
	throw std::runtime_error("not implemented yet");
}

double gcContent(const std::string& kmer) {
	throw std::runtime_error("not implemented yet");
}

std::string kmerAfterError(const std::string& kmer, size_t pos, ErrorType type) {
	throw std::runtime_error("not implemented yet");
}

std::unordered_map<size_t, size_t> countReadLengths(const std::vector<std::string>& readFiles) {
	throw std::runtime_error("not implemented yet");
}

Dataset::Dataset(GenomeType genomeType, size_t genomeSize, const std::vector<std::string>& readFiles) :
		_genomeType(genomeType), _genomeSize(genomeSize), _readFiles(readFiles) {
	_readLengths = countReadLengths(readFiles);
}

ReferenceDataset::ReferenceDataset(GenomeType genomeType, size_t genomeSize, const std::vector<std::string>& readFiles,
		const std::string& referenceGenomePath) :
		Dataset(genomeType, genomeSize, readFiles), _referenceGenomePath(referenceGenomePath) {
	_referenceGenome = io::readReferenceGenome(referenceGenomePath);
}

} // end of namespace seq_correct::util
} // end of namespace seq_correct
