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

#include <stddef.h>
#include <unordered_map>
#include "../../external/cereal/types/unordered_map.hpp"
#include "../../seq_correct.hpp"

namespace seq_correct {
namespace profile {

class MotifTreeEntry {
public:
	MotifTreeEntry();
	void reset();

	std::unordered_map<ErrorType, size_t> numErrors;
	std::unordered_map<ErrorType, double> zScore;

	template<class Archive>
	void serialize(Archive & archive) {
		archive(numErrors, zScore); // serialize things by passing them to the archive
	}
};

} // end of namespace seq_correct::profile
} // end of namespace seq_correct
