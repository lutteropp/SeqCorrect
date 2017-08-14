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

/**
 * @brief This header includes all other SeqCorrect headers (*.hpp).
 * This makes it easy to integrate the code as a library, as all
 * symbols of SeqCorrect are available after including this header.
 */

#pragma once

#include "coverage/coverage_bias.hpp"
#include "correction/error_correction.hpp"
#include "counting/fm_count.hpp"
#include "eval/error_evaluation_data.hpp"
#include "profile/error_profile.hpp"
#include "pusm/pusm.hpp"
#include "io/sequence_io.hpp"
#include "util/util.hpp"
#include "util/dataset.hpp"
#include "eval/evaluation.hpp"
#include "eval/metrics.hpp"
#include "io/bam_iterator.hpp"
#include "io/read_with_alignments.hpp"
#include "external/bloom_filter.hpp"
#include "kmer/classification.hpp"
#include "util/enums.hpp"
