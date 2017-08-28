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
#include "../../external/cereal/types/unordered_map.hpp"
#include "../../external/cereal/types/vector.hpp"

#include "../error_profile.hpp"

#include "motif_tree.hpp"

namespace seq_correct {
namespace profile {

static const int MAX_MOTIF_SIZE = 6;


class MotifErrorProfile : public ErrorProfileUnit {
public:
	MotifErrorProfile(Matcher &kmerCounter);
	virtual std::unordered_map<ErrorType, double> getErrorProbabilities(const Read &read, size_t positionInRead);
	virtual std::unordered_map<ErrorType, double> getKmerErrorProbabilities(const std::string &kmer, size_t positionInKmer);
	virtual void loadErrorProfile(const std::string &filepath, Matcher &counter);
	virtual void storeErrorProfile(const std::string &filepath);
	virtual void plotErrorProfile();

	virtual void reset();
	virtual void check(const std::vector<Correction>& corrections, const Read &originalRead);

	virtual void finalize();

	double findMostSignificantZScore(const ErrorType &type, const std::string &sequence, int posInSequence);

	template<class Archive>
	void serialize(Archive & archive) {
		archive(motifTree, finalized); // serialize things by passing them to the archive
	}
protected:
	virtual std::unordered_map<ErrorType, double> getErrorProbabilitiesFinalized(const Read &read, size_t positionInRead);
	virtual std::unordered_map<ErrorType, double> getErrorProbabilitiesFinalized(const std::string &kmer, size_t positionInKmer);
private:
	void updateMotifData(int positionInSequence, const std::string &sequence, const ErrorType &type);
	void computeZScores();

	MotifTree motifTree;
	bool finalized;
	Matcher &counter;
};

} // end of namespace seq_correct::profile
} // end of namespace seq_correct
