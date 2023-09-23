#include <cstdint>
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

#ifndef COUNT_CU_SKETCH_H
#define COUNT_CU_SKETCH_H

class CountCuSketch {

private:
	size_t w, d;
	std::vector<std::shared_ptr<std::vector<int32_t>>> matrix;
	std::vector<int> seeds;
	std::vector<int> sign_seeds;

public:
	CountCuSketch(size_t w, size_t d): w(w), d(d)
	{
		srand(10);
		for (size_t i = 0; i < d; ++i)
		{
			matrix.push_back (std::make_shared <std::vector<int32_t>> (1 << w, 0));
			seeds.push_back(rand());
		}
	}

	~CountCuSketch(){};

	int64_t InsertSketch(packet item)
	{
		size_t p = item.hash (seeds[0]) & ((1 << w) - 1);
		int32_t minval = matrix[0]->at(p);

		for (size_t i = 1; i < d; ++i)
		{
			size_t p = item.hash (seeds[i]) & ((1 << w) - 1);
			if (matrix[i]->at(p) < minval) minval = matrix[i]->at(p);
		}

		for (size_t i = 0; i < d; ++i)
		{
			size_t p = item.hash (seeds[i]) & ((1 << w) - 1);
			int s = (matrix[i]->at(p) == minval) ? 1 : 0;
			int64_t value = matrix[i]->at(p) + s;
			matrix[i]->at(p) = value;
		}

		return minval + 1;
	}

	int64_t QuerySketch(packet item)
	{
		size_t p = item.hash (seeds[0]) & ((1 << w) - 1);
		int32_t minval = matrix[0]->at(p);

		for (size_t i = 1; i < d; ++i)
		{
			size_t p = item.hash (seeds[i]) & ((1 << w) - 1);
			if (matrix[i]->at(p) < minval) minval = matrix[i]->at(p);
		}
		return minval;
	}
};

#endif // COUNT_CU_SKETCH_H

