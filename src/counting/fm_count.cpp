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

bool isSelfReverseComplement(const std::string& kmer) {
	if (kmer.size() % 2 == 1) {
		return false;
	}
	bool same = true;
	for (size_t i = 0; i < kmer.size() / 2; ++i) {
		if (kmer[i] != util::reverseComplementChar(kmer[kmer.size() - i - 1])) {
			same = false;
			break;
		}
	}
	return same;
}

uint16_t FMIndexMatcher::countKmer(const std::string& kmer) {
	if (useBuffer && buffer.find(kmer) != buffer.end()) {
		return buffer[kmer];
	}
	size_t count = sdsl::count(fmIndex, kmer.begin(), kmer.end());
	if (isSelfReverseComplement(kmer)) {
		count /= 2;
	}
	if (useBuffer) {
		buffer[kmer] = count;
	}
	return count;
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

NaiveBufferedMatcher::NaiveBufferedMatcher(const std::string& filename, size_t k, bool revCompExtra) {
	io::ReadInput input;

	std::string filenameNaive = filename + "." + std::to_string(k) + ".naive";
	std::ifstream test(filenameNaive);
	if (test.good()) {
		std::string kmer;
		uint16_t count;
		while (test >> kmer >> count) {
			buffer[kmer] = count;
		}
	} else {
		input.openFile(filename);
		double minProgress = 0;
		while (input.hasNext()) {
			if (input.progress() > minProgress) {
				std::cout << input.progress() << "\n";
				minProgress += 1;
			}
			io::Read read = input.readNext(true, false, false);
			for (size_t i = 0; i < read.seq.size() - k; ++i) {
				std::string kmer = read.seq.substr(i, k);
				if (buffer.find(kmer) != buffer.end()) {
					buffer[kmer]++;
				} else {
					buffer[kmer] = 1;
				}

				if (revCompExtra) {
					kmer = util::reverseComplementString(kmer);
					if (buffer.find(kmer) != buffer.end()) {
						buffer[kmer]++;
					} else {
						buffer[kmer] = 1;
					}
				}
			}
		}

		std::ofstream testOut(filenameNaive);
		for (auto kv : buffer) {
			testOut << kv.first << " " << kv.second << "\n";
		}
		testOut.close();
	}
}

uint16_t NaiveBufferedMatcher::countKmer(const std::string& kmer) {
	uint16_t count = buffer[kmer];
	if (isSelfReverseComplement(kmer)) {
		count /= 2;
	}
	return count;
}

} // end of namespace seq_correct::counting
} // end of namespace seq_correct
