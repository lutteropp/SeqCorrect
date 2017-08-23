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
#include <unordered_set>
#include <iostream>

#include "fm_count.hpp"
#include "../util/util.hpp"

namespace seq_correct {
namespace counting {

void createReadsOnly(const std::string& readsPath) {
	// create the readsOnly file if it doesn't exist yet
	std::ifstream testInput(readsPath + ".readsOnly.txt");
	if (testInput.good()) {
		testInput.close();
		return; // file already exists, do nothing.
	}
	std::cout << "Creating " << readsPath + ".readsOnly.txt" << "...\n";
	io::ReadInput reader;
	reader.openFile(readsPath);
	std::ofstream writer(readsPath + ".readsOnly.txt");
	while (reader.hasNext()) {
		std::string seq = reader.readNext(true, false, false).seq;
		writer << seq << "$" << util::reverseComplementString(seq) << "$"; // use '$' as a delimiter
	}
	writer.close();
	std::cout << "Successfully created " << readsPath + ".readsOnly.txt" << "...\n";
}

uint16_t FMIndexMatcher::countKmer(const std::string& kmer) {
	if (useBuffer && buffer.find(kmer) != buffer.end()) {
		return buffer[kmer];
	}
	size_t count = sdsl::count(fmIndex, kmer.begin(), kmer.end());
	if (util::isSelfReverseComplement(kmer)) {
		count /= 2;
	}
	if (useBuffer) {
		buffer[kmer] = count;
	}
	return count;
}

FMIndexMatcher::FMIndexMatcher(const std::string& readsPath) {
	createReadsOnly(readsPath);
	std::string filename = readsPath + ".readsOnly.txt";

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
	useBuffer = false;
}

void FMIndexMatcher::enableBuffer() {
	useBuffer = true;
}

void FMIndexMatcher::disableBuffer() {
	useBuffer = false;
}

FMIndexMatcherMulti::FMIndexMatcherMulti(const std::string& filename) :
		FMIndexMatcher(filename) {
	std::ifstream infile(filename);
	std::istreambuf_iterator<char> it(infile);
	size_t actId = 0;
	while (it != std::istreambuf_iterator<char>()) {
		if (*it == '$') {
			actId++;
		}
		documentID.push_back(actId / 2); // divided by 2 because the reverse-complement it put afterwards
		++it;
	}
}

std::vector<uint16_t> FMIndexMatcherMulti::countKmerMultiPositions(const std::string& kmer,
		bool returnEmptyIfDoubleOccs) {
	auto locations = sdsl::locate(fmIndex, kmer.begin(), kmer.end());
	std::sort(locations.begin(), locations.end());
	std::vector<uint16_t> res;
	std::unordered_set<uint16_t> documents;
	for (size_t i = 0; i < locations.size(); ++i) {
		if (returnEmptyIfDoubleOccs && documents.find(documentID[locations[i]]) != documents.end()) {
			return std::vector<uint16_t>();
		}
		if (returnEmptyIfDoubleOccs) {
			documents.insert(documentID[locations[i]]);
		}
		res.push_back(documentID[locations[i]]);
	}
	return res;
}

} // end of namespace seq_correct::counting
} // end of namespace seq_correct
