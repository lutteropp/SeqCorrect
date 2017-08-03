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
#include "external/tclap/CmdLine.h"
#include "util/genome_type.hpp"

using namespace seq_correct;
using namespace util;

void createReadsOnly(const std::string& pathToOriginalReads) {
	// create the readsOnly file if it doesn't exist yet
	std::ifstream testInput(pathToOriginalReads + ".readsOnly.txt");
	if (testInput.good()) {
		testInput.close();
		return; // file already exists, do nothing.
	}
	io::ReadInput reader;
	reader.openFile(pathToOriginalReads);
	std::ofstream writer(pathToOriginalReads + ".readsOnly.txt");
	while (reader.hasNext()) {
		writer << reader.readNext(true, false, false).seq << "\n";
	}
	writer.close();
}

void cmd_demo(const std::string& outputPath) {
	// Specify an example Ebola Illumina dataset
	util::Dataset dataset(util::GenomeType::LINEAR, 18959,
			"/home/sarah/Documents/Master Thesis Topic Extension/thesis_zipped/SOFTWARE_AND_DATA/data/Simulated Datasets/Ebola/Illumina/ebola_illumina_simulated.fastq",
			"/home/sarah/Documents/Master Thesis Topic Extension/thesis_zipped/SOFTWARE_AND_DATA/data/Simulated Datasets/Ebola/Illumina/ebola_illumina_simulated.fastq.readsOnly.txt");

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
	counting::FMIndexMatcher matcher(dataset.getReadFilepath());
	size_t count = matcher.countKmerNoRC("ACGGT");
	std::cout << "The kmer 'ACGGT' occurs " << count << " times in the read dataset\n";
	count = matcher.countKmer("ACGGT");
	std::cout << "The kmer 'ACGGT' or its reverse-complement occurs " << count << " times in the read dataset\n";

	// compute the coverage biases
	pusm::PerfectUniformSequencingModel pusm(dataset.getGenomeType(), dataset.getGenomeSize(),
			dataset.getReadLengths());
	coverage::CoverageBiasUnitSingle biasUnit;
	biasUnit.preprocess(21, dataset.getReadFilepath(), matcher, pusm);
	biasUnit.printMedianCoverageBiases();
}

void cmd_correct(const std::string& pathToOriginalReads, GenomeType genomeType, const std::string& outputPath, size_t genomeSize) {
	createReadsOnly(pathToOriginalReads);
	std::string pathToReadsOnly = pathToOriginalReads + ".readsOnly.txt";
	counting::FMIndexMatcher fm(pathToReadsOnly);
	std::unordered_map<size_t, size_t> readLengths = countReadLengths(pathToOriginalReads);
	pusm::PerfectUniformSequencingModel pusm(genomeType, genomeSize, readLengths);
	coverage::CoverageBiasUnitMulti biasUnit;
	// TODO: compute coverage biases
	// TODO: correct the reads
	// TODO: print the result
}

void cmd_eval(const std::string& pathToOriginalReads, const std::string& pathToCorrectedReads,
		const std::string& pathToGenome, const std::string& outputPath) {

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

	try {
		TCLAP::CmdLine cmd("SeqCorrect - an error correction toolkit for next-generation whole-genome sequencing reads",
				' ', "0.1");
		TCLAP::ValueArg<std::string> readsArg("r", "reads", "Path to the sequencing reads in FASTA or FASTQ format",
				true, "", "string");
		cmd.add(readsArg);
		TCLAP::ValueArg<std::string> inputArg("i", "input", "Path to the corrected reads in FASTA or FASTQ format",
				false, "", "string");
		cmd.add(inputArg);
		TCLAP::ValueArg<std::string> genomeArg("g", "genome", "Path to the reference genome in FASTA format", false, "",
				"string");
		cmd.add(genomeArg);
		TCLAP::ValueArg<size_t> genomeSizeArg("s", "size", "Estimated genome size", false, 0, "unsigned int");
		cmd.add(genomeSizeArg);
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
	} catch (TCLAP::ArgException &e) // catch any exceptions
	{
		std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
		return 1;
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
	if (correctMode && !(circular || linear)) {
		throw std::runtime_error("Either --circular or --linear is required");
	}
	if (!correctMode && (circular || linear)) {
		throw std::runtime_error("Genome type only required for --correct mode!");
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

	if (demoMode) {
		cmd_demo(outputPath);
	} else if (correctMode) {
		if (circular) {
			genomeType = GenomeType::CIRCULAR;
		} else {
			genomeType = GenomeType::LINEAR;
		}
		cmd_correct(pathToOriginalReads, genomeType, outputPath, genomeSize);
	} else if (evalMode) {
		cmd_eval(pathToOriginalReads, pathToCorrectedReads, pathToGenome, outputPath);
	}
}
