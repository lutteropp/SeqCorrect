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
#include <condition_variable>
#include <functional>
#include <mutex>
#include <queue>
#include <vector>
#include <iostream>
#include <cassert>
#include <thread>

namespace seq_correct {
namespace util {

template<class T> class ProducerConsumerPattern {
public:
	ProducerConsumerPattern(size_t maxDataQueueSize, std::function<double(std::vector<T>&, size_t)> produce, std::function<void(std::vector<T>&, size_t)> consume);
	void run(size_t numProducers, size_t numConsumers);
private:
	void produce(size_t producerId);
	void consume(size_t consumerId, size_t myProducer);

	std::function<double(std::vector<T>&, size_t)> prepareData;
	std::function<void(std::vector<T>&, size_t)> processData;

	std::vector<std::mutex> mtx;
	std::vector<std::condition_variable> cvConsume;
	std::vector<std::condition_variable> cvProduce;
	std::vector<std::queue<std::vector<T> > > dataQueue;

	std::vector<bool> done;

	size_t _maxDataQueueSize;
};


template<class T> ProducerConsumerPattern<T>::ProducerConsumerPattern(size_t maxDataQueueSize,
		std::function<double(std::vector<T>&, size_t)> produce,
		std::function<void(std::vector<T>&, size_t)> consume) {
	_maxDataQueueSize = maxDataQueueSize;
	prepareData = produce;
	processData = consume;
}

// assumes that the called functions are threadsafe
template<class T> void ProducerConsumerPattern<T>::run(size_t numProducers, size_t numConsumers) {
	std::vector<std::thread> threads;

	cvConsume = std::vector<std::condition_variable>(numProducers);
	cvProduce = std::vector<std::condition_variable>(numProducers);
	mtx = std::vector<std::mutex>(numProducers);
	dataQueue = std::vector<std::queue<std::vector<T> > >(numProducers);
	done = std::vector<bool>(numProducers, false);

	assert(numConsumers % numProducers == 0);
	size_t consumersPerProducer = numConsumers / numProducers;

	for (size_t i = 0; i < numProducers; ++i) {
		threads.push_back(std::thread(&ProducerConsumerPattern<T>::produce, this, i));
	}
	for (size_t i = 0; i < numConsumers; ++i) {
		threads.push_back(std::thread(&ProducerConsumerPattern<T>::consume, this, i, i / consumersPerProducer));
	}
	for (size_t i = 0; i < threads.size(); ++i) {
		threads[i].join();
	}
	std::cout << "Threads finished.\n";
}

template<class T> void ProducerConsumerPattern<T>::produce(size_t producerId) {
	size_t actReads = 0;
	double minProgress = 1;
	while (!done[producerId]) {
		std::vector<T> buffer;
		double progress = prepareData(buffer, producerId);
		actReads += buffer.size();
		std::unique_lock<std::mutex> lck(mtx[producerId]);

		cvProduce[producerId].wait(lck, [this, &producerId] {return dataQueue[producerId].size() <= _maxDataQueueSize;});
		dataQueue[producerId].push(buffer);
		cvConsume[producerId].notify_one();

		if (progress >= minProgress) {
			std::cout << "Producer " + std::to_string(producerId) + ": " + std::to_string(progress) + " \%\n";
			minProgress += 1;
		}
		if (progress == 100) {
			done[producerId] = true;
			cvConsume[producerId].notify_all();
		}
	}
	std::cout << "Producer " + std::to_string(producerId) + ": finished producing.\n";
}

template<class T> void ProducerConsumerPattern<T>::consume(size_t consumerId, size_t myProducer) {
	while (!done[myProducer] || !dataQueue.empty()) {
		std::unique_lock<std::mutex> lck(mtx[myProducer]);
		cvConsume[myProducer].wait(lck, [this, &myProducer] {return (!dataQueue[myProducer].empty()) || done[myProducer];});
		{
			if (dataQueue[myProducer].empty()) {
				break;
			}
			std::vector<T> buffer = dataQueue[myProducer].front();
			dataQueue[myProducer].pop();
			lck.unlock();
			cvProduce[myProducer].notify_one();
			processData(buffer, consumerId);
		}
	}
	std::cout << "Consumer " + std::to_string(consumerId) + ": finished consuming.\n";
}

} // end of namespace seq_correct::util
} // end of namespace seq_correct
