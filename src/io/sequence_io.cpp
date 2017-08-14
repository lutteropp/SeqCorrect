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
#include <experimental/filesystem>
#include <sstream>
#include <iostream>
#include <ios>
#include "sequence_io.hpp"

namespace seq_correct {
namespace io {

/**
 * Read a genome in FASTA format. Convert all its bases to upper case.
 * @param filepath Path to the FASTA file containig the reference genome
 */
std::string readReferenceGenome(const std::string& filepath) {
	if (filepath.empty()) {
		throw std::runtime_error("Genome filepath is empty!");
	}
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

/**
 * Check if the file contains more reads.
 */
bool ReadInput::hasNext() {
	int len = file.tellg(); // current position in the ifstream
	std::string line;
	getline(file, line);
	bool res = file.good();
	file.seekg(len, std::ios_base::beg); // return back to position before reading the line
	return res;
}

/**
 * Open a FASTA or FASTQ file containing the reads.
 * @param filepath The file containing the reads
 */
void ReadInput::openFile(const std::string& filepath) {
	file.open(filepath);

	if (!file.good()) {
		throw std::runtime_error("Could not open file: " + filepath);
	}

	// count number of reads
	numReadsTotal = 0;
	std::ifstream temp(filepath);
	std::string line;
	while (std::getline(temp, line)) {
		if (line[0] == '@' || line[0] == '>') {
			numReadsTotal++;
		}
	}
	readsLeft = numReadsTotal;
	temp.close();
}

double ReadInput::progress() {
	return ((numReadsTotal - readsLeft) / ((double) numReadsTotal)) * 100;
}

/**
 * Read the next read in the input file.
 * @param readSequence Read the DNA sequence
 * @param readQuality Read the quality scores
 * @param readName Read the sequence name
 */
Read ReadInput::readNext(bool readSequence, bool readQuality, bool readName) {
	bool isFASTA = true;
	std::string line;
	std::getline(file, line);
	if (line[0] == '@') {
		isFASTA = false;
	} else if (line[0] != '>') {
		throw std::runtime_error("the sequence is neither in FASTA nor in FASTQ format");
	}
	if (readQuality && isFASTA) {
		throw std::runtime_error("cannot read quality scores");
	}
	Read read;
	if (readName) {
		read.name = line.substr(1, std::string::npos);
	}
	std::getline(file, line);
	if (readSequence) {
		read.seq = line;
		std::transform(read.seq.begin(), read.seq.end(), read.seq.begin(), ::toupper);
	}
	if (!isFASTA) {
		std::getline(file, line); // the '+'
		std::getline(file, line); // the quality scores
		if (readQuality) {
			read.qual = line;
		}
	}
	readsLeft--;
	return read;
}

/**
 * Create a new file to write the reads to it.
 * @param filepath The path to the newly created file
 */
bool ReadOutput::createFile(const std::string& filepath) {
	if (std::experimental::filesystem::exists(filepath)) {
		throw std::runtime_error("the file does already exist");
	}
	file.open(filepath);
	return file.good();
}

/**
 * Write a read to the output file.
 * @param read The read
 */
void ReadOutput::write(const Read& read) {
	if (!file.good()) {
		throw std::runtime_error("can not access output file - did you create one?");
	}
	std::stringstream buffer;
	if (read.qual.empty()) { // FASTA format
		buffer << ">" << read.name << "\n" << read.seq << "\n";
	} else { // FASTQ format
		buffer << "@" << read.name << "\n" << read.seq << "\n+\n" << read.qual << "\n";
	}
	file << buffer.str();
}

void ReadOutput::close() {
	file.close();
}

} // end of namespace seq_correct::io
} // end of namespace seq_correct
