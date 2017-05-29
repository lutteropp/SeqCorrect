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

namespace io {

/**
 * Reads a reference genome from a file.
 * @param filepath path to a FASTA file containing the reference genome
 */
std::string readReferenceGenome(const std::string& filepath);

class Read {
public:
	Read() :
			_seq(""), _qual(""), _name("") {
	}
	Read(std::string&& seq);
	void setSeq(const std::string& seq);
	void setQual(const std::string& qual);
	void setName(const std::string& name);
private:
	std::string _seq;
	std::string _qual;
	std::string _name;
};

// TODO: Don't forget to convert all DNA sequences to upper case
class ReadInput {
public:
	bool openFile(const std::string& filepath);
	Read readNext(bool readSequence = true, bool readQuality = false, bool readName = false);
	bool hasNext();
private:
	std::ifstream _file;
};

class ReadOutput {
public:
	bool openFile(const std::string& filepath);
	bool write(const Read& read);
private:
	std::ofstream _file;
};

} // end of namespace seq_correct::io
