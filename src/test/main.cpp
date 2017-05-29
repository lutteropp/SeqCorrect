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

#include <iostream>
#include <string>
#include <gtest/gtest.h>
#include "../seq_correct.hpp"

using namespace helper;

TEST(HelperTest, ReverseComplementString) {
    ASSERT_EQ("A", reverseComplementString("T"));
    ASSERT_EQ("C", reverseComplementString("G"));
    ASSERT_EQ("G", reverseComplementString("C"));
    ASSERT_EQ("T", reverseComplementString("A"));
    ASSERT_EQ("N", reverseComplementString("N"));
    ASSERT_EQ("AT", reverseComplementString("AT")); // reverse-complement of itself
    ASSERT_EQ("CG", reverseComplementString("CG")); // reverse-complement of itself
    ASSERT_EQ("GA", reverseComplementString("TC"));
    ASSERT_EQ("TG", reverseComplementString("CA"));
}

int main(int argc, char **argv) {
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
