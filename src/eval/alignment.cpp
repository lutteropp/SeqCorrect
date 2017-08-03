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

#include <bitset>
#include <stdexcept>
#include <functional>
#include <utility>
#include <limits>
#include <iostream>
#include "alignment.hpp"

/*
 * TODO: Idea: In case of multiple optimal alignments... maybe return all of them and then later take the one that
 * has the smallest distance to the corresponding ground truth read?
 */

namespace seq_correct {
namespace eval {

size_t coordinate(size_t m, size_t i, size_t j) {
	return i * (m + 1) + j;
}

size_t min3(const size_t& a, const size_t& b, const size_t& c) {
	return std::min(a, std::min(b, c));
}

/*
 * In case of multiple optimal alignments, only return one of them.
 */
void backtrack(size_t i, size_t j, const std::string& originalRead, const std::string& correctedRead,
		const std::vector<size_t>& M, const std::vector<std::bitset<3> >& info,
		std::vector<std::pair<size_t, util::ErrorType> >& res) {
	if (i == 0 && j == 0) { // check if backtracking has ended
		return;
	}
	size_t n = originalRead.size();
	size_t m = correctedRead.size();
	std::function<size_t(size_t, size_t)> coord = std::bind(coordinate, m, std::placeholders::_1,
			std::placeholders::_2);
	size_t actCoord = coord(i, j);
	if (info[actCoord].count() > 1) {
		std::cout << "WARNING: These sequences have multiple optimal alignments... \n" << "\tOriginal Read: "
				<< originalRead << "\n\tCorrectedRead: " << correctedRead << "\n";
	}
	if (info[actCoord][0] == 1) { // match or mismatch
		if (originalRead[i] != correctedRead[j]) {
			res.push_back(std::make_pair(i, util::inferSubstitutionFrom(correctedRead[j])));
		}
		backtrack(i-1, j-1, originalRead, correctedRead, M, info, res);
	} else if (info[actCoord][1] == 1) { // insertion into original read, i.e. we had a deletion within the original
		res.push_back(std::make_pair(i, util::inferDeletionOf(correctedRead[j])));
		backtrack(i, j-1, originalRead, correctedRead, M, info, res);
	} else if (info[actCoord][2] == 1) { // deletion into original read, i.e. we had an insertion within the original
		res.push_back(std::make_pair(i, util::ErrorType::INSERTION));
		backtrack(i-1, j, originalRead, correctedRead, M, info, res);
	} else {
		throw std::runtime_error("This should not happen.");
	}

}

void smoothenMultidels(std::vector<std::pair<size_t, util::ErrorType> >& res) {
	size_t i = 0;
	if (res.size() < 2) return;
	while (i < res.size() - 1) {
		if (res[i].first == res[i+1].first) {
			if (util::isGapErrorType(res[i].second) && util::isGapErrorType(res[i+1].second)) {
				res[i].second = util::ErrorType::MULTIDEL;
				res.erase(res.begin() + i+1);
				--i;
			}
		}
		++i;
	}
}

std::vector<std::pair<size_t, util::ErrorType> > align(const std::string& originalRead,
		const std::string& correctedRead) {
	size_t mis = 1; // mismatch penalty
	size_t ins = 1; // insertion into original read
	size_t del = 1; // deletion from original read
	std::vector<std::pair<size_t, util::ErrorType> > res;

	size_t n = originalRead.size();
	size_t m = correctedRead.size();
	size_t size = (n + 1) * (m + 1);
	std::function<size_t(size_t, size_t)> coord = std::bind(coordinate, m, std::placeholders::_1,
			std::placeholders::_2);
	std::vector<size_t> M(size);
	std::vector<std::bitset<3> > info(size);

	// initialization
	M[coord(0, 0)] = 0;
	for (size_t i = 1; i < n + 1; ++i) {
		M[coord(i, 0)] = std::numeric_limits<size_t>::infinity();
	}
	for (size_t j = 1; j < m + 1; ++j) {
		M[coord(0, j)] = std::numeric_limits<size_t>::infinity();
	}

	// matrix fillup
	for (size_t i = 1; i < n + 1; ++i) {
		for (size_t j = 1; j < m + 1; ++j) {
			// check matrix M
			size_t mismatch = (originalRead[i - 1] != correctedRead[i - 1]) ? mis : 0;
			size_t min = min3(M[coord(i - 1, j - 1)] + mismatch, M[coord(i, j - 1)] + ins, M[coord(i - 1, j) + del]);
			if (M[coord(i - 1, j - 1)] + mismatch == min) {
				info[coord(i, j)][0] = 1;
			}
			if (M[coord(i, j - 1)] + ins == min) {
				info[coord(i, j)][1] = 1;
			}
			if (M[coord(i - 1, j)] + del == min) {
				info[coord(i, j)][2] = 1;
			}
			M[coord(i, j)] = min;
		}
	}

	backtrack(n, m, originalRead, correctedRead, M, info, res);
	smoothenMultidels(res);
	return res;
}

} // end of namespace eval
} // end of namespace seq_correct
