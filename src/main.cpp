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

int main() {
	std::vector<int> vec = { 0, 1, 2, 3, 4, 5, 6 };
	auto itBegin = vec.begin();
	auto itEnd = vec.end();
	auto itRBegin = vec.rbegin();
	auto itREnd = vec.rend();
	std::cout << "*itBegin: " << *itBegin << std::endl;
	std::cout << "*itEnd: " << *itEnd << std::endl;
	std::cout << "*itRBegin: " << *itRBegin << std::endl;
	std::cout << "*itREnd: " << *itREnd << std::endl;

	std::cout << "*itBegin + 3: " << *(itBegin + 3) << std::endl;
	std::cout << "*itRBegin + 3: " << *(itRBegin + 3) << std::endl;
	std::cout << "*itREnd + 3: " << *(itREnd + 3) << std::endl;

	std::cout << "*itEnd - 3: " << *(itEnd - 3) << std::endl;
	std::cout << "*itRBegin - 3: " << *(itRBegin - 3) << std::endl;
	std::cout << "*itREnd - 2: " << *(itREnd - 2) << std::endl;
}
