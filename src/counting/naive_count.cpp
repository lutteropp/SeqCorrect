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

#include "naive_count.hpp"

namespace seq_correct {
namespace counting {

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
	if (util::isSelfReverseComplement(kmer)) {
		count /= 2;
	}
	return count;
}

} // end of namespace seq_correct::counting
} // end of namespace seq_correct
