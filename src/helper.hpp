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

#include "sequence_io.hpp"

namespace helper {

	enum class ErrorType {
		INSERTION, SUB_OF_A, SUB_OF_C, SUB_OF_G, SUB_OF_T, DEL_OF_A, DEL_OF_C, DEL_OF_G, DEL_OF_T, MULTIDEL
	};

	std::string reverseComplementString(const std::string& text);
	io::Read reverseComplementRead(const io::Read& read);
	double gcContent(const std::string& kmer);
	std::string kmerAfterError(const std::string& kmer, size_t pos, ErrorType type);

} // end of namespace helper
