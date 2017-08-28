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
#include <cassert>
#include <fstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "../external/cereal/archives/binary.hpp"
#include "../seq_correct.hpp"

namespace seq_correct {
namespace profile {

using namespace util;
using namespace io;
using namespace eval;
using namespace counting;

class ErrorProfileUnit {
public:
	virtual ~ErrorProfileUnit() {};

	virtual std::vector<std::unordered_map<ErrorType, double> > getReadErrorProbabilities(const Read &read) {
		std::vector<std::unordered_map<ErrorType, double> > res;
		for (size_t i = 0; i < read.seq.size(); ++i) {
			res.push_back(getErrorProbabilities(read, i));
		}
		return res;
	}

	virtual std::vector<std::unordered_map<ErrorType, double> > getKmerErrorProbabilities(const std::string &kmer) {
			std::vector<std::unordered_map<ErrorType, double> > res;
			for (size_t i = 0; i < kmer.size(); ++i) {
				res.push_back(getKmerErrorProbabilities(kmer, i));
			}
			return res;
		}

	virtual std::vector<std::unordered_map<ErrorType, double> > getReadErrorProbabilitiesPartial(const Read &read, size_t from, size_t to) {
			std::vector<std::unordered_map<ErrorType, double> > res;
			assert(from < read.sequence.size());
			assert(to < read.sequence.size());
			for (size_t i = from; i <= to; ++i) {
				res.push_back(getErrorProbabilities(read, i));
			}
			return res;
		}

	virtual void learnErrorProfileFromFiles(const std::string &correctionsFile, double acceptProb = 1.0) {
		reset();
		std::ifstream infile;
		infile.open(correctionsFile, std::ios::binary);
		if (!infile.good()) {
			throw std::runtime_error("The file " + correctionsFile + " does not exist!");
		}
		size_t minProgress = 1;
		size_t n;
		Read cr;
		cereal::BinaryInputArchive iarchive(infile);
		iarchive(n);
		for (size_t i = 0; i < n; ++i) {
			iarchive(cr);
			check(cr, acceptProb);
			double progress = (double) i * 100 / n;
			if (progress >= minProgress) {
				std::cout << progress << "%\n";
				minProgress++;
			}
		}
		infile.close();
		finalize();
	}

	virtual void learnErrorProfileFromFilesAligned(const std::string &correctionsFile, double acceptProb = 1.0) {
		reset();
		std::ifstream infile;
		infile.open(correctionsFile, std::ios::binary);
		if (!infile.good()) {
			throw std::runtime_error("The file " + correctionsFile + " does not exist!");
		}
		size_t minProgress = 0;
		unsigned long long n;
		ReadWithAlignments cra;
		cereal::BinaryInputArchive iarchive(infile);
		iarchive(n);
		for (unsigned long long i = 0; i < n; ++i) {
			iarchive(cra);
			checkAligned(cra, acceptProb);
			double progress = (double) i * 100 / n;
			if (progress >= minProgress) {
				std::cout << progress << "%\n";
				minProgress++;
			}
		}
		infile.close();
		finalize();
	}

	virtual std::unordered_map<ErrorType, double> getErrorProbabilities(const Read &read, size_t positionInRead) = 0;
	virtual std::unordered_map<ErrorType, double> getKmerErrorProbabilities(const std::string &kmer, size_t positionInKmer) = 0;
	virtual void loadErrorProfile(const std::string &filepath, Matcher &counter) = 0;
	virtual void storeErrorProfile(const std::string &filepath) = 0;
	virtual void plotErrorProfile() = 0;

	virtual void reset() = 0;
	virtual void check(const Read &corrRead, double acceptProb = 1.0) = 0;
	virtual void checkAligned(const ReadWithAlignments &corrRead, double acceptProb = 1.0) = 0;
	virtual void finalize() = 0;
protected:
	virtual std::unordered_map<ErrorType, double> getErrorProbabilitiesFinalized(const Read &read, size_t positionInRead) = 0;
	virtual std::unordered_map<ErrorType, double> getErrorProbabilitiesFinalized(const std::string &kmer, size_t positionInKmer) = 0;
};

} // end of namespace seq_correct::profile
} // end of namespace seq_correct
