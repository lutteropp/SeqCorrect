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

#include "kmer_evaluation.hpp"
#include "../kmer/classification.hpp"
#include "../kmer/hash_classifier.hpp"
#include "../io/bam_iterator.hpp"
#include "../counting/counting.hpp"

namespace seq_correct {
namespace eval {

using namespace classification;

void classifyKmersVariants(size_t k, GenomeType genomeType, const std::string& pathToOriginalReads,
		const std::string& genomeFilepath, counting::Matcher& fmReads, counting::Matcher& fmGenome);

std::vector<KmerType> computeTrueTypes(size_t k, const std::string& sequence, counting::Matcher& fmGenome) {
	std::vector<KmerType> trueTypes;
	for (size_t i = 0; i < sequence.size() - k; ++i) {
		std::string kmer = sequence.substr(i, k);
		size_t trueCount = fmGenome.countKmer(kmer);
		KmerType trueType;
		if (trueCount == 0) {
			trueType = KmerType::UNTRUSTED;
		} else if (trueCount == 1) {
			trueType = KmerType::UNIQUE;
		} else {
			trueType = KmerType::REPEAT;
		}
		trueTypes.push_back(trueType);
	}
	return trueTypes;
}

std::vector<KmerType> computeTrueTypes(size_t k, const std::string& sequence, HashClassifierGenome& clsfyGenome) {
	std::vector<KmerType> trueTypes;
	for (size_t i = 0; i < sequence.size() - k; ++i) {
		std::string kmer = sequence.substr(i, k);
		KmerType trueType = clsfyGenome.classifyKmer(kmer);
		trueTypes.push_back(trueType);
	}
	return trueTypes;
}

double computeMedianKmerCount(const std::vector<uint16_t>& counts) {
	double medianCount = 0;
	std::vector<uint16_t> kmerCountsSorted;
	for (size_t i = 0; i < counts.size(); ++i) {
		kmerCountsSorted.push_back(counts[i]);
	}
	std::sort(kmerCountsSorted.begin(), kmerCountsSorted.end());
	if (kmerCountsSorted.size() % 2 == 1) {
		medianCount = (double) kmerCountsSorted[(kmerCountsSorted.size() + 1) / 2];
	} else {
		medianCount = ((double) kmerCountsSorted[kmerCountsSorted.size() / 2]
				+ (double) kmerCountsSorted[kmerCountsSorted.size() / 2 + 1]) / 2;
	}
	if (medianCount < 3) {
		medianCount = 3;
	}
	return medianCount;
}

std::vector<uint16_t> computeKmerCountsRead(size_t k, const std::string& sequence, counting::Matcher& fmReads) {
	std::vector<uint16_t> kmerCounts;
	for (size_t i = 0; i < sequence.size() - k; ++i) {
		std::string kmer = sequence.substr(i, k);
		uint16_t count = fmReads.countKmer(kmer);
		kmerCounts.push_back(count);
	}
	return kmerCounts;
}

double computeAverageBias(size_t k, const std::string& sequence, coverage::CoverageBiasUnitMulti& biasUnit,
		const std::string& pathToOriginalReads, counting::Matcher& fmReads, pusm::PerfectUniformSequencingModel& pusm) {
	double biasSum = 0;
	size_t biasNum = 0;
	for (size_t i = 0; i < sequence.size() - k; ++i) {
		std::string kmer = sequence.substr(i, k);
		biasSum += biasUnit.computeCoverageBias(kmer, pathToOriginalReads, fmReads, pusm);
		biasNum++;
	}
	return biasSum / biasNum;
}

void update(const std::vector<KmerType>& trueTypes, const std::vector<KmerType>& predictedTypes,
		KmerEvaluationData& data, size_t& numCorrect, size_t& numWrong) {
	bool correctlyClassified = true;
	for (size_t i = 0; i < trueTypes.size(); ++i) {
		KmerType trueType = trueTypes[i];
		KmerType predictedType = predictedTypes[i];
		if (trueType != predictedType) {
			correctlyClassified = false;
		}
		data.update(trueType, predictedType);
	}
	if (correctlyClassified) {
		numCorrect++;
	} else {
		/*std::cout << "Read K-mer-decomposition:\n";
		 for (size_t i = 0; i < trueTypes.size(); ++i) {
		 std::cout << "  " << kmerTypeToString(trueTypes[i]) << " -> " << kmerTypeToString(predictedTypes[i])
		 << ": " << kmerCounts[i] << "\n";
		 }*/

		numWrong++;
	}
}

void printKmerEvaluationData(const eval::KmerEvaluationData& evalData) {
	for (KmerType type : KmerTypeIterator()) {
		std::cout << kmerTypeToString(type) << ":\n";
		std::cout << "  TP:          " << eval::truePositives(type, evalData) << "\n";
		std::cout << "  TN:          " << eval::trueNegatives(type, evalData) << "\n";
		std::cout << "  FP:          " << eval::falsePositives(type, evalData) << "\n";
		std::cout << "  FN:          " << eval::falseNegatives(type, evalData) << "\n";
		std::cout << "  Accuracy:    " << eval::computeAccuracy(type, evalData) << "\n";
		std::cout << "  Precision:   " << eval::computePrecision(type, evalData) << "\n";
		std::cout << "  Recall:      " << eval::computeRecall(type, evalData) << "\n";
		std::cout << "  Specificity: " << eval::computeSpecificity(type, evalData) << "\n";
		std::cout << "  Sensitivity: " << eval::computeSensitivity(type, evalData) << "\n";
		std::cout << "  Gain:        " << eval::computeGain(type, evalData) << "\n";
		std::cout << "  F1-score:    " << eval::computeF1Score(type, evalData) << "\n";
	}
	std::cout << "\n";
	std::cout << "Unweighted Average F1-score: " << eval::computeUnweightedAverageF1Score(evalData) << "\n";
	std::cout << "NMI-score:                   " << eval::computeNMIScore(evalData) << "\n";

	std::cout << "\n Confusion Matrix:\n";
	for (KmerType e1 : KmerTypeIterator()) {
		for (KmerType e2 : KmerTypeIterator()) {
			double percentage = (double) evalData.getEntry(e1, e2) / evalData.sumTruth(e1);
			std::cout << "[" << util::kmerTypeToString(e1) << "][" << util::kmerTypeToString(e2) << "]: "
					<< evalData.getEntry(e1, e2) << " = " << percentage * 100 << "%" << "\n";
		}
	}
}

void classifyKmersVariants(size_t k, GenomeType genomeType, const std::string& pathToOriginalReads,
		const std::string& genomeFilepath, counting::Matcher& fmReads, counting::Matcher& fmGenome) {
	std::string genome = io::readReferenceGenome(genomeFilepath);
	std::unordered_map<size_t, size_t> readLengths = countReadLengths(pathToOriginalReads);
	pusm::PerfectUniformSequencingModel pusm(genomeType, genome.size(), readLengths);
	coverage::CoverageBiasUnitMulti biasUnit;

	size_t numCorrectThesis = 0;
	size_t numWrongThesis = 0;
	size_t numCorrectThesisRead1 = 0;
	size_t numWrongThesisRead1 = 0;
	size_t numCorrectThesisRead2 = 0;
	size_t numWrongThesisRead2 = 0;
	size_t numCorrectRead = 0;
	size_t numWrongRead = 0;

	KmerEvaluationData dataThesis;
	KmerEvaluationData dataThesisRead1;
	KmerEvaluationData dataThesisRead2;
	KmerEvaluationData dataRead;

	io::ReadInput readInput(pathToOriginalReads);
	double expectedCount = pusm.expectedCount(k).expectation;

	double minProgress = 0.0;
	while (readInput.hasNext()) {
		std::string seq = readInput.readNext(true, false, false).seq;

		std::vector<uint16_t> kmerCounts = computeKmerCountsRead(k, seq, fmReads);
		double medianCount = computeMedianKmerCount(kmerCounts);
		std::vector<KmerType> trueTypes = computeTrueTypes(k, seq, fmGenome);

		double biasReadAverage = computeAverageBias(k, seq, biasUnit, pathToOriginalReads, fmReads, pusm);
		double biasReadK = biasUnit.computeCoverageBias(k, (double) util::countGC(seq) / seq.size(),
				pathToOriginalReads, fmReads, pusm);

		std::vector<KmerType> predictedTypesThesis(seq.size() - k);
		std::vector<KmerType> predictedTypesThesisRead1(seq.size() - k);
		std::vector<KmerType> predictedTypesThesisRead2(seq.size() - k);
		std::vector<KmerType> predictedTypesRead(seq.size() - k);
#pragma omp parallel for
		for (size_t i = 0; i < seq.size() - k; ++i) {
			std::string kmer = seq.substr(i, k);

			size_t observedCount = kmerCounts[i];
			double biasKmer = biasUnit.computeCoverageBias(kmer, pathToOriginalReads, fmReads, pusm);

			predictedTypesThesis[i] = classifyKmer(observedCount, expectedCount, biasKmer);
			predictedTypesThesisRead1[i] = classifyKmer(observedCount, expectedCount, biasReadAverage);
			predictedTypesThesisRead2[i] = classifyKmer(observedCount, expectedCount, biasReadK);
			predictedTypesRead[i] = classifyKmerReadBased(k, i, kmerCounts, medianCount, seq);
		}

		update(trueTypes, predictedTypesThesis, dataThesis, numCorrectThesis, numWrongThesis);
		update(trueTypes, predictedTypesThesisRead1, dataThesisRead1, numCorrectThesisRead1, numWrongThesisRead1);
		update(trueTypes, predictedTypesThesisRead2, dataThesisRead2, numCorrectThesisRead2, numWrongThesisRead2);
		update(trueTypes, predictedTypesRead, dataRead, numCorrectRead, numWrongRead);

		if (readInput.progress() > minProgress) {
			minProgress += 1;
			std::cout << readInput.progress() << "\n";
		}
	}

	std::cout << "numCorrectThesis: " << numCorrectThesis << " = "
			<< (double) numCorrectThesis * 100 / (numCorrectThesis + numWrongThesis) << "% \n";
	std::cout << "numWrongThesis: " << numWrongThesis << " = "
			<< (double) numWrongThesis * 100 / (numCorrectThesis + numWrongThesis) << "% \n";

	std::cout << "numCorrectThesisRead1: " << numCorrectThesisRead1 << " = "
			<< (double) numCorrectThesisRead1 * 100 / (numCorrectThesisRead1 + numWrongThesisRead1) << "% \n";
	std::cout << "numWrongThesisRead1: " << numWrongThesisRead1 << " = "
			<< (double) numWrongThesisRead1 * 100 / (numCorrectThesisRead1 + numWrongThesisRead1) << "% \n";

	std::cout << "numCorrectThesisRead2: " << numCorrectThesisRead2 << " = "
			<< (double) numCorrectThesisRead2 * 100 / (numCorrectThesisRead2 + numWrongThesisRead2) << "% \n";
	std::cout << "numWrongThesisRead2: " << numWrongThesisRead2 << " = "
			<< (double) numWrongThesisRead2 * 100 / (numCorrectThesisRead2 + numWrongThesisRead2) << "% \n";

	std::cout << "numCorrectRead: " << numCorrectRead << " = "
			<< (double) numCorrectRead * 100 / (numCorrectRead + numWrongRead) << "% \n";
	std::cout << "numWrongRead: " << numWrongRead << " = "
			<< (double) numWrongRead * 100 / (numCorrectRead + numWrongRead) << "% \n";

	std::cout << "Variant Thesis k-mer classification: \n";
	printKmerEvaluationData(dataThesis);

	std::cout << "\nVariant ThesisRead1 k-mer classification: \n";
	printKmerEvaluationData(dataThesisRead1);

	std::cout << "\nVariant ThesisRead2 k-mer classification: \n";
	printKmerEvaluationData(dataThesisRead2);

	std::cout << "\nVariant Read k-mer classification: \n";
	printKmerEvaluationData(dataRead);
}

void classifyKmers(size_t k, GenomeType genomeType, const std::string& pathToOriginalReads,
		const std::string& genomeFilepath, counting::Matcher& fmReads, counting::Matcher& fmGenome) {
	std::string genome = io::readReferenceGenome(genomeFilepath);
	std::unordered_map<size_t, size_t> readLengths = countReadLengths(pathToOriginalReads);
	pusm::PerfectUniformSequencingModel pusm(genomeType, genome.size(), readLengths);
	coverage::CoverageBiasUnitMulti biasUnit;
	biasUnit.preprocess(k, pathToOriginalReads, fmReads, pusm);
	HashClassifier clsfyReads(fmReads, pusm, biasUnit, pathToOriginalReads);
	HashClassifierGenome clsfyGenome(fmGenome);

	size_t numCorrectThesis = 0;
	size_t numWrongThesis = 0;

	KmerEvaluationData dataThesis;

	io::ReadInput readInput(pathToOriginalReads);

	double minProgress = 0.0;
	while (readInput.hasNext()) {
		std::string seq = readInput.readNext(true, false, false).seq;

		std::vector<KmerType> trueTypes = computeTrueTypes(k, seq, clsfyGenome);

		std::vector<KmerType> predictedTypesThesis(seq.size() - k);
#pragma omp parallel for
		for (size_t i = 0; i < seq.size() - k; ++i) {
			std::string kmer = seq.substr(i, k);
			predictedTypesThesis[i] = clsfyReads.classifyKmer(kmer);
		}

		update(trueTypes, predictedTypesThesis, dataThesis, numCorrectThesis, numWrongThesis);

		if (readInput.progress() > minProgress) {
			minProgress += 1;
			std::cout << readInput.progress() << "\n";
		}
	}

	std::cout << "numCorrectThesis: " << numCorrectThesis << " = "
			<< (double) numCorrectThesis * 100 / (numCorrectThesis + numWrongThesis) << "% \n";
	std::cout << "numWrongThesis: " << numWrongThesis << " = "
			<< (double) numWrongThesis * 100 / (numCorrectThesis + numWrongThesis) << "% \n";

	std::cout << "Variant Thesis k-mer classification: \n";
	printKmerEvaluationData(dataThesis);
}

void eval_kmers(size_t k, GenomeType genomeType, const std::string& pathToOriginalReads,
		const std::string& pathToGenome) {
	counting::NaiveBufferedMatcher fmReads(pathToOriginalReads, k, true);
	//counting::NaiveBufferedMatcher fmGenome(pathToGenome, k, true);

	//counting::FMIndexMatcher fmReads(pathToOriginalReads);
	counting::FMIndexMatcher fmGenome(pathToGenome);

	std::cout << "k-mer size used: " << k << "\n";
	eval::classifyKmers(k, genomeType, pathToOriginalReads, pathToGenome, fmReads, fmGenome);

	//eval::classifyKmersVariants(k, genomeType, pathToOriginalReads, pathToGenome, fmReads, fmGenome);
}

} // end of namespace seq_correct::eval
} // end of namespace seq_correct
