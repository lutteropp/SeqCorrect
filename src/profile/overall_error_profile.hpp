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
#include <map>
#include <string>
#include <unordered_map>
#include <utility>
#include "../external/cereal/types/unordered_map.hpp"
#include "../external/cereal/types/utility.hpp"

#include "error_profile.hpp"

namespace seq_correct {
namespace profile {

using namespace util;
using namespace io;

class OverallErrorProfile : public ErrorProfileUnit {
public:
	OverallErrorProfile();
	virtual std::unordered_map<ErrorType, double> getErrorProbabilities(const Read &read, size_t positionInRead);
	virtual std::unordered_map<ErrorType, double> getKmerErrorProbabilities(const std::string &kmer, size_t positionInKmer);
	virtual void loadErrorProfile(const std::string &filepath, Matcher &counter);
	virtual void storeErrorProfile(const std::string &filepath);
	virtual void plotErrorProfile();

	virtual void reset();
	virtual void check(const Read &corrRead, double acceptProb = 1.0);
	virtual void checkAligned(const ReadWithAlignments &corrRead, double acceptProb = 1.0);

	virtual void finalize();

	double getOverallErrorRateCurrentBase();
	double getOverallErrorRateNextGap();

	template<class Archive>
	void serialize(Archive & archive) {
		archive(counts, substitutionMatrix, totalCount, noncorrectBases, deletedBases); // serialize things by passing them to the archive
	}
protected:
	virtual std::unordered_map<ErrorType, double> getErrorProbabilitiesFinalized(const Read &read, size_t positionInRead);
	virtual std::unordered_map<ErrorType, double> getErrorProbabilitiesFinalized(const std::string &kmer, size_t positionInKmer);
private:
	std::unordered_map<ErrorType, size_t> counts;
	std::unordered_map<ErrorType, double> counts_finalized;
	std::unordered_map<std::pair<char, char>, size_t, pairhash> substitutionMatrix;
	std::unordered_map<std::pair<char, char>, double, pairhash> substitutionMatrix_finalized;
	size_t totalCount;
	size_t noncorrectBases;
	size_t deletedBases;

	double overallErrorProbCurrent;
	double overallErrorProbNext;

	bool finalized;
};

} // end of namespace seq_correct::profile
}
 // end of namespace seq_correct
