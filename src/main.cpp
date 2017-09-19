/*
 * main.cpp
 *
 *  Created on: May 30, 2017
 *      Author: sarah
 */

#include "seq_correct.hpp"
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <fstream>
#include <cstdlib>
#include "external/tclap/CmdLine.h"

using namespace seq_correct::util;
using namespace seq_correct;

void cmd_demo(const std::string& outputPath) {
	// Specify an example Ebola Illumina dataset
	util::Dataset dataset(util::GenomeType::LINEAR, 18959,
			"/home/sarah/Documents/Master Thesis Topic Extension/thesis_zipped/SOFTWARE_AND_DATA/data/Simulated Datasets/Ebola/Illumina/ebola_illumina_simulated.fastq");

	// Read the reads and write them
	/*io::ReadInput reader;
	 reader.openFile(dataset.getReadFilepath());
	 io::ReadOutput writer;
	 writer.createFile("temp.txt");
	 while (reader.hasNext()) {
	 io::Read read = reader.readNext(true, false, true);
	 writer.write(read);
	 }*/

	// count the k-mer "ACGGT" in all the reads
	counting::FMIndexMatcher matcher(dataset.getPathToOriginalReads());
	size_t count = matcher.countKmer("ACGGT");
	std::cout << "The kmer 'ACGGT' or its reverse-complement occurs " << count << " times in the read dataset\n";

	// compute the coverage biases
	pusm::PerfectUniformSequencingModel pusm(dataset.getGenomeType(), dataset.getGenomeSize(),
			dataset.getReadLengths());
	coverage::CoverageBiasUnitSingle biasUnit;
	biasUnit.preprocess(21, dataset.getPathToOriginalReads(), matcher, pusm);
	biasUnit.printMedianCoverageBiases();
}

void cmd_correct(size_t k, const std::string& pathToOriginalReads, GenomeType genomeType, const std::string& outputPath,
		size_t genomeSize, correction::CorrectionAlgorithm algo) {
	counting::FMIndexMatcher fm(pathToOriginalReads);
	std::unordered_map<size_t, size_t> readLengths = countReadLengths(pathToOriginalReads);
	pusm::PerfectUniformSequencingModel pusm(genomeType, genomeSize, readLengths);
	coverage::CoverageBiasUnitMulti biasUnit;
	correction::correctReads(pathToOriginalReads, algo, fm, pusm, biasUnit, outputPath, k);
}

void cmd_eval(size_t k, GenomeType genomeType, const std::string& pathToOriginalReads,
		const std::string& pathToCorrectedReads, const std::string& pathToGenome, const std::string& outputPath) {
	std::cout << pathToCorrectedReads << "\n\n";
	eval::eval_corrections(k, genomeType, pathToOriginalReads, pathToCorrectedReads, pathToGenome, outputPath);
	//eval::eval_kmers(k, genomeType, pathToOriginalReads, pathToGenome);
}

int main(int argc, char* argv[]) {
	std::string pathToOriginalReads;
	GenomeType genomeType;
	bool circular = false;
	bool linear = false;
	bool demoMode; // load an example dataset, good for quick testing and code demonstrations
	bool correctMode;
	bool evalMode;
	size_t genomeSize = 0;
	std::string pathToCorrectedReads;
	std::string pathToGenome;
	std::string outputPath;
	std::string algoName;
	correction::CorrectionAlgorithm algo = correction::CorrectionAlgorithm::SIMPLE_KMER;
	size_t k;

	try {
		TCLAP::CmdLine cmd("SeqCorrect - an error correction toolkit for next-generation whole-genome sequencing reads",
				' ', "0.1");
		TCLAP::ValueArg<std::string> readsArg("r", "reads", "Path to the sequencing reads in FASTA or FASTQ format",
				false, "", "string");
		cmd.add(readsArg);
		TCLAP::ValueArg<std::string> inputArg("i", "input", "Path to the corrected reads in FASTA or FASTQ format",
				false, "", "string");
		cmd.add(inputArg);
		TCLAP::ValueArg<std::string> genomeArg("g", "genome", "Path to the reference genome in FASTA format", false, "",
				"string");
		cmd.add(genomeArg);
		TCLAP::ValueArg<size_t> genomeSizeArg("s", "size", "Estimated genome size", false, 0, "unsigned int");
		cmd.add(genomeSizeArg);
		TCLAP::ValueArg<size_t> kmerSizeArg("k", "kmer", "K-mer size to use", false, 17, "unsigned int");
		cmd.add(kmerSizeArg);
		TCLAP::ValueArg<std::string> outputArg("o", "output", "Path to the output file", true, "", "string");
		cmd.add(outputArg);
		TCLAP::SwitchArg circularArg("c", "circular", "Circular genome", false);
		cmd.add(circularArg);
		TCLAP::SwitchArg linearArg("l", "linear", "Linear genome", false);
		cmd.add(linearArg);
		TCLAP::SwitchArg demoArg("d", "demo", "Demonstration mode", false);
		cmd.add(demoArg);
		TCLAP::SwitchArg evalArg("e", "eval", "Evaluation mode", false);
		cmd.add(evalArg);
		TCLAP::SwitchArg correctArg("n", "normal", "Error Correction mode", false);
		cmd.add(correctArg);

		TCLAP::ValueArg<std::string> algoArg("a", "algo", "Error correction algorithm", false, "", "string");
		cmd.add(algoArg);

		cmd.parse(argc, argv);

		pathToOriginalReads = readsArg.getValue();
		pathToCorrectedReads = inputArg.getValue();
		pathToGenome = genomeArg.getValue();
		outputPath = outputArg.getValue();
		circular = circularArg.getValue();
		linear = linearArg.getValue();
		demoMode = demoArg.getValue();
		correctMode = correctArg.getValue();
		evalMode = evalArg.getValue();
		genomeSize = genomeSizeArg.getValue();
		algoName = algoArg.getValue();
		k = kmerSizeArg.getValue();
	} catch (TCLAP::ArgException &e) // catch any exceptions
	{
		std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
		return 1;
	}

	if (correctMode && algoName.empty()) {
		throw std::runtime_error(
				"You must specify an algorithm! Accepted choices are: NONE, SIMPLE_KMER, ADAPTIVE_KMER, FULL_MSA, PARTIAL_MSA, SUFFIX_TREE");
	}

	if (correctMode && !algoName.empty()) {
		if (algoName == "NONE") {
			algo = correction::CorrectionAlgorithm::NONE;
		} else if (algoName == "SIMPLE_KMER") {
			algo = correction::CorrectionAlgorithm::SIMPLE_KMER;
		} else if (algoName == "ADAPTIVE_KMER") {
			algo = correction::CorrectionAlgorithm::ADAPTIVE_KMER;
		} else if (algoName == "FULL_MSA") {
			algo = correction::CorrectionAlgorithm::FULL_MSA;
		} else if (algoName == "PARTIAL_MSA") {
			algo = correction::CorrectionAlgorithm::PARTIAL_MSA;
		} else if (algoName == "SUFFIX_TREE") {
			algo = correction::CorrectionAlgorithm::SUFFIX_TREE;
		} else {
			throw std::runtime_error(
					"Unknown algorithm! Accepted choices are: NONE, SIMPLE_KMER, ADAPTIVE_KMER, FULL_MSA, PARTIAL_MSA, SUFFIX_TREE");
		}
	}

	if (correctMode && (genomeSize == 0)) {
		throw std::runtime_error("You must specify an estimated genome size > 0!");
	}

	if (!correctMode && (genomeSize != 0)) {
		throw std::runtime_error("--size only allowed in --normal mode");
	}

	if (demoMode + correctMode + evalMode != 1) {
		throw std::runtime_error("You must specify exactly one out of --demo, --normal, and --eval");
	}
	if (circular && linear) {
		throw std::runtime_error("A genome cannot be both circular and linear");
	}
	if ((correctMode || evalMode) && !(circular || linear)) {
		throw std::runtime_error("Either --circular or --linear is required");
	}
	if (!(correctMode || evalMode) && (circular || linear)) {
		throw std::runtime_error("Genome type only required for --correct and --eval mode!");
	}
	if (evalMode && pathToGenome.empty()) {
		throw std::runtime_error("Path to genome is missing");
	}
	if (evalMode && pathToCorrectedReads.empty()) {
		throw std::runtime_error("Path to corrected reads is missing");
	}
	if (!evalMode && (!pathToCorrectedReads.empty())) {
		throw std::runtime_error("Path to corrected reads only required for --eval option");
	}
	if (!evalMode && (!pathToGenome.empty())) {
		throw std::runtime_error("Path to genome only required for --eval option");
	}
	if (demoMode
			&& (!pathToCorrectedReads.empty() || !pathToGenome.empty() || !pathToOriginalReads.empty() || circular
					|| linear)) {
		throw std::runtime_error(
				"Please do not specify any other parameters than --output, when using the --demo mode");
	}

	if (circular) {
		genomeType = GenomeType::CIRCULAR;
	} else {
		genomeType = GenomeType::LINEAR;
	}

	if (demoMode) {
		cmd_demo(outputPath);
	} else if (correctMode) {
		cmd_correct(k, pathToOriginalReads, genomeType, outputPath, genomeSize, algo);
	} else if (evalMode) {
		cmd_eval(k, genomeType, pathToOriginalReads, pathToCorrectedReads, pathToGenome, outputPath);
	}
}
