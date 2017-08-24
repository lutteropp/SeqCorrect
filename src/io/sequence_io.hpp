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

#include <string>
#include <vector>
#include <fstream>

namespace seq_correct {
namespace io {

/**
 * Reads a reference genome from a file.
 * @param filepath path to a FASTA file containing the reference genome
 */
std::string readReferenceGenome(const std::string& filepath);

struct Read {
public:
	Read() :
			seq(""), qual(""), name("") {
	}
	Read(std::string&& seq);

	std::string seq;
	std::string qual;
	std::string name;
};

class ReadInput {
public:
	ReadInput(const std::string& filepath);
	Read readNext(bool readSequence = true, bool readQuality = false, bool readName = true);
	bool hasNext();
	double progress();
private:
	void openFile(const std::string& filepath);
	std::ifstream file;
	size_t numReadsTotal;
	size_t readsLeft;
};

class ReadOutput {
public:
	bool createFile(const std::string& filepath);
	void write(const Read& read);
	void close();
private:
	std::ofstream file;
};

} // end of namespace seq_correct::io
} // end of namespace seq_correct
