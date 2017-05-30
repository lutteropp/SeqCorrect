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

#include "pusm.hpp"
#include "info.hpp"
#include <stdexcept>

namespace seq_correct{
namespace pusm {

PusmData expectedCountLinear(size_t genomeSize, const std::unordered_map<size_t, size_t>& readLengths, size_t k) {
	throw std::runtime_error("not implemented yet");
}

PusmData expectedCountCircular(size_t genomeSize, const std::unordered_map<size_t, size_t>& readLengths, size_t k) {
	throw std::runtime_error("not implemented yet");
}

PusmData PerfectUniformSequencingModel::expectedCount(size_t k) {
	if (_pusmBuffer.find(k) != _pusmBuffer.end()) {
		return _pusmBuffer[k];
	}

	PusmData res;
	if (_type == info::GenomeType::CIRCULAR) {
		res = expectedCountCircular(_genomeSize, _readLengths, k);
	} else {
		res = expectedCountLinear(_genomeSize, _readLengths, k);
	}
	_pusmBuffer[k] = res;
	return res;
}

} // end of namespace seq_correct::pusm
} // end of namespace seq_correct
