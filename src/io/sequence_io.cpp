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
#include <algorithm>
#include "sequence_io.hpp"

namespace seq_correct {
namespace io {

/**
 * Read a genome in FASTA format. Convert all its bases to upper case.
 * @param filepath Path to the FASTA file containig the reference genome
 */
std::string readReferenceGenome(const std::string& filepath) {
	std::string genome;
	std::ifstream infile(filepath);
	std::string line;
	std::getline(infile, line);
	if (line.empty() || line[0] != '>') {
		throw std::runtime_error("not a FASTA file");
	}
	while (std::getline(infile, line)) {
		if (!line.empty() && line[0] == '>') {
			throw std::runtime_error("multiple sequences in this file");
		}
		genome += line;
	}
	std::transform(genome.begin(), genome.end(), genome.begin(), ::toupper);
	return genome;
}

bool ReadInput::hasNext() {
	throw std::runtime_error("not implemented yet");
}

bool ReadInput::openFile(const std::string& filepath) {
	throw std::runtime_error("not implemented yet");
}

// TODO: Don't forget to convert all DNA sequences to upper case
Read ReadInput::readNext(bool readSequence, bool readQuality, bool readName) {
	throw std::runtime_error("not implemented yet");
}

bool ReadOutput::openFile(const std::string& filepath) {
	throw std::runtime_error("not implemented yet");
}

bool ReadOutput::write(const Read& read) {
	throw std::runtime_error("not implemented yet");
}

} // end of namespace seq_correct::io
} // end of namespace seq_correct