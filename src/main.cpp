/*
 * main.cpp
 *
 *  Created on: May 30, 2017
 *      Author: sarah
 */

#include "seq_correct.hpp"
#include <vector>
#include <algorithm>
#include <iostream>

using namespace seq_correct;

int main() {
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
	pusm::PerfectUniformSequencingModel pusm(dataset.getGenomeType(), dataset.getGenomeSize(), dataset.getReadLengths());
	coverage::CoverageBiasUnit biasUnit;
	biasUnit.preprocess(21, dataset.getReadFilepath(), matcher, pusm);
	biasUnit.printMedianCoverageBiases();
}
