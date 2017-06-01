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

#include <stdexcept>

#include "fm_count.hpp"
#include "../util/util.hpp"

namespace seq_correct {
namespace counting {

size_t FMIndexMatcher::countKmer(const std::string& kmer) {
	return countKmerNoRC(kmer) + countKmerNoRC(util::reverseComplementString(kmer));
}

size_t FMIndexMatcher::countKmer(const external::ConstStringPtr& kmer) {
	return countKmer(kmer.toString());
}

size_t FMIndexMatcher::countKmerNoRC(const std::string& kmer) {
	return sdsl::count(fmIndex, kmer.begin(), kmer.end());
}

size_t FMIndexMatcher::countKmerNoRC(const external::ConstStringPtr& kmer) {
	return countKmerNoRC(kmer.toString());
}

FMIndexMatcher::FMIndexMatcher(const std::string& filename) {
	std::ifstream infile(filename);
	if (!infile.good()) {
		throw std::runtime_error("This file does not exist! " + filename);
	}

	std::string index_suffix = ".fm9";
	std::string index_file = filename + index_suffix;
	if (!load_from_file(fmIndex, index_file)) {
		std::cout << "No FM index found. Constructing FM index..." << std::endl;
		construct(fmIndex, filename, 1); // generate index
		store_to_file(fmIndex, index_file); // save it
		std::cout << "Index construction complete, index requires " << size_in_mega_bytes(fmIndex) << " MiB."
				<< std::endl;
	}
}

/*FMIndexMatcher buildIndex(const std::string &text) {
 throw std::runtime_error("not implemented yet");
 }

 FMIndexMatcher buildIndex(const std::vector<std::string> &filepaths) {
 throw std::runtime_error("not implemented yet");
 }

 FMIndexMatcher buildIndex(const io::Read &read) {
 throw std::runtime_error("not implemented yet");
 }

 FMIndexMatcher buildIndex(const std::vector<io::Read> &reads) {
 throw std::runtime_error("not implemented yet");
 }

 FMIndexMatcher loadIndex(const std::string &filepath) {
 throw std::runtime_error("not implemented yet");
 }

 bool storeIndex(const std::string &filepath, const FMIndexMatcher &index) {
 throw std::runtime_error("not implemented yet");
 }*/

} // end of namespace seq_correct::counting
} // end of namespace seq_correct
