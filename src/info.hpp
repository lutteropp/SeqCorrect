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

#include <unordered_map>
#include <vector>
#include <string>

#include "sequence_io.hpp"

namespace info {

enum class GenomeType {
	CIRCULAR, LINEAR
};

class Dataset {
public:
	Dataset(GenomeType genomeType, size_t genomeSize, const std::vector<std::string>& readFiles);
private:
	GenomeType _genomeType;
	size_t _genomeSize;
	std::unordered_map<size_t, size_t> _readLengths;
	std::vector<std::string> &_readFiles;
};

class ReferenceDataset: public Dataset {
public:
	ReferenceDataset(GenomeType genomeType, size_t genomeSize, const std::vector<std::string>& readFiles,
			const std::string& referenceGenomePath);
private:
	std::string _referenceGenomePath;
	std::string _referenceGenome;
};

} // end of namespace info
