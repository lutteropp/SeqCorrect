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

#include "dataset.hpp"
#include "../io/sequence_io.hpp"

namespace seq_correct {
namespace util {

std::unordered_map<size_t, size_t> countReadLengths(const std::string& readFilepath) {
	std::unordered_map<size_t, size_t> res;
	io::ReadInput reader;
	reader.openFile(readFilepath);
	while (reader.hasNext()) {
		io::Read read = reader.readNext(true, false, false);
		if (res.find(read.seq.size()) == res.end()) {
			res[read.seq.size()] = 1;
		} else {
			res[read.seq.size()]++;
		}
	}
	return res;
}

Dataset::Dataset(GenomeType genomeType, size_t genomeSize, const std::string& readFilepath,
		const std::string& readsOnlyFilepath) :
		genomeType(genomeType), genomeSize(genomeSize), readFilepath(readFilepath), readsOnlyFilepath(readsOnlyFilepath) {
	readLengths = countReadLengths(readFilepath);
}

size_t Dataset::getGenomeSize() {
	return genomeSize;
}

GenomeType Dataset::getGenomeType() {
	return genomeType;
}

std::string Dataset::getReadFilepath() {
	return readFilepath;
}

std::unordered_map<size_t, size_t> Dataset::getReadLengths() {
	return readLengths;
}

std::string Dataset::getReadsOnlyFilepath() {
	return readsOnlyFilepath;
}

ReferenceDataset::ReferenceDataset(GenomeType genomeType, size_t genomeSize, const std::string& readFilepath,
		const std::string& referenceGenomePath, const std::string& readsOnlyFilepath) :
		Dataset(genomeType, genomeSize, readFilepath, readsOnlyFilepath), referenceGenomePath(referenceGenomePath) {
	referenceGenome = io::readReferenceGenome(referenceGenomePath);
}

std::string ReferenceDataset::getReferenceGenome() {
	return referenceGenome;
}

std::string ReferenceDataset::getReferenceGenomePath() {
	return referenceGenomePath;
}

} // end of namespace seq_correct::util
} // end of namespace seq_cirrect
