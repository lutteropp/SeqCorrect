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
std::string readReferenceGenome(const std::string &filepath);

// TODO: This has so much code duplication... how to do it the proper C++ 11 way? I smell templates, but I don't know how...

// TODO: Don't forget to convert all DNA sequences to upper case

class FASTARead {
public:
	FASTARead(std::string &&seq);
private:
	std::string _seq;
};

class FASTQRead : public FASTARead {
public:
	FASTQRead(std::string &&seq, std::string &&qual);
	FASTQRead(FASTARead &&read, std::string &&qual);
private:
	std::string _qual;
};

class NamedFASTARead : public FASTARead {
public:
	NamedFASTARead(std::string &&seq, std::string &&name);
	NamedFASTARead(FASTARead &&read, std::string &&name);
private:
	std::string _name;
};

class NamedFASTQRead : public FASTQRead {
public:
	NamedFASTQRead(std::string &&seq, std::string &&qual, std::string &&name);
	NamedFASTQRead(FASTQRead &&read, std::string &&name);
private:
	std::string _name;
};

class FASTAReader {
public:
	bool openFile(const std::string &filepath);
	FASTARead readNext();
	bool hasNext();
private:
	std::ifstream _file;
};

class FASTQReader {
public:
	bool openFile(const std::string &filepath);
	FASTQRead readNext();
	bool hasNext();
private:
	std::ifstream _file;
};

class NamedFASTAReader {
public:
	bool openFile(const std::string &filepath);
	NamedFASTARead readNext();
	bool hasNext();
private:
	std::ifstream _file;
};

class NamedFASTQReader {
public:
	bool openFile(const std::string &filepath);
	NamedFASTQRead readNext();
	bool hasNext();
private:
	std::ifstream _file;
};

class FASTAWriter {
public:
	bool openFile(const std::string &filepath);
	bool write(const FASTARead &read);
private:
	std::ifstream _file;
};

class FASTQWriter {
public:
	bool openFile(const std::string &filepath);
	bool write(const FASTQRead &read);
private:
	std::ifstream _file;
};

class NamedFASTAWriter {
public:
	bool openFile(const std::string &filepath);
	bool write(const NamedFASTARead &read);
private:
	std::ifstream _file;
};

class NamedFASTQWriter {
public:
	bool openFile(const std::string &filepath);
	bool write(const NamedFASTQRead &read);
private:
	std::ifstream _file;
};

} // end of namespace seq_correct::io
