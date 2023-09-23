#include <iostream>
#include <map>
#include <vector>
#include <unordered_map>
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

#ifndef COUNTMIN_PERFECT_H
#define COUNTMIN_PERFECT_H

class CountPerfect {

private:
	std::unordered_map<packet, int32_t> memory;

public:
	CountPerfect(size_t w, size_t d)
	{
	}

	~CountPerfect(){};

	int64_t InsertSketch(packet item)
	{
		if (memory.find(item) != memory.end())
			memory[item] ++;
		else
			memory[item] = 1;

		return memory[item];
	}

	int64_t QuerySketch(packet item)
	{
		return memory[item];
	}
};

#endif // COUNTMIN_PERFECT_H


