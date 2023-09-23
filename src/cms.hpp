#include <iostream>
#include <map>
#include <vector>
#include <cmath>
#include <algorithm>
#include <queue>
#include <utility>
#include <string>
#include <iomanip>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <limits>
#include <cstddef>
#include <memory>
#include "packet.hpp"

using namespace std;

#ifndef COUNTMIN_SKETCH_H
#define COUNTMIN_SKETCH_H

class CountminSketch {

private:
	size_t w, d;
	std::vector<std::shared_ptr<std::vector<int32_t>>> matrix;
	std::vector<int> seeds;

public:
	CountminSketch(size_t w, size_t d): w(w), d(d)
	{
		srand(10);
		for (size_t i = 0; i < d; ++i)
		{
			matrix.push_back (std::make_shared <std::vector<int32_t>> (1 << w, 0));
			seeds.push_back(rand());
		}
	}

	~CountminSketch(){};

	int64_t InsertSketch(packet item)
	{
		int64_t minval = std::numeric_limits<int>::max();
		for (size_t i = 0; i < d; ++i)
		{
			size_t p = item.hash (seeds[i]) & ((1 << w) - 1);
			matrix[i]->at(p) = matrix[i]->at(p) + 1;
			if (matrix[i]->at(p) < minval) minval = matrix[i]->at(p);
		}
		return minval;
	}

	int64_t QuerySketch(packet item)
	{
		int64_t minval = std::numeric_limits<long>::max();

		for (size_t i = 0; i < d; ++i)
		{
			size_t p = item.hash (seeds[i]) & ((1 << w) - 1);
			if (matrix[i]->at(p) < minval) minval = matrix[i]->at(p);
		}
		return minval;
	}
};

#endif // COUNTMIN_SKETCH_H

