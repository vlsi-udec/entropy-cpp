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

#ifndef COUNT_SKETCH_H
#define COUNT_SKETCH_H

class CountSketch {

private:
	size_t w, d;
	std::vector<std::shared_ptr<std::vector<int32_t>>> matrix;
	std::vector<int> seeds;
	std::vector<int> sign_seeds;

public:
	CountSketch(size_t w, size_t d): w(w), d(d)
	{
		srand(10);
		for (size_t i = 0; i < d; ++i)
		{
			matrix.push_back (std::make_shared <std::vector<int32_t>> (1 << w, 0));
			seeds.push_back(rand());
			sign_seeds.push_back(rand());
		}
	}

	~CountSketch(){};

	int64_t InsertSketch (packet item)
	{
		int64_t values[this->d];
		for (size_t i = 0; i < d; ++i)
		{
			size_t p = item.hash (seeds[i]) & ((1 << w) - 1);
			int s = (item.hash (sign_seeds[i]) & 0x01) ? -1 : 1;
			int64_t value = matrix[i]->at(p) + s;
			matrix[i]->at(p) = value;
			values[i] = s * value;
		}

		nth_element(values, values + d / 2, values + d);
		return values[d / 2];
	}

	int64_t QuerySketch(packet item)
	{
		int64_t values[this->d];
		for (size_t i = 0; i < d; ++i)
		{
			size_t p = item.hash (seeds[i]) & ((1 << w) - 1);
			int s = (item.hash (sign_seeds[i]) & 0x01) ? -1 : 1;
			values[i] = s * matrix[i]->at(p);
		}
		// return the median (4.3.2 The median trick, "ESTIMATING THE NUMBER OF DISTINCT ELEMENTS", page 18)
		nth_element(values, values + d / 2, values + d);
		return values[d / 2];
	}
};

#endif // COUNT_SKETCH_H

