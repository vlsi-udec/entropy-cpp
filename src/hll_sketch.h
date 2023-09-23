#include <iostream>
#include <errno.h>
#include <string>
#include <math.h>
#include <stdlib.h>
#include <map>
#include <vector>
#include <string.h>
#include <algorithm>
#include <unordered_set>
#include <bits/stdc++.h>
#include "parametros.h"
#include "hllpp.hpp"

#ifndef HLL_SKETCH_H_INCLUDED
#define HLL_SKETCH_H_INCLUDED

using namespace std;

inline uint count_zeros (uint64_t value, uint precision)
{
	return (__builtin_clzll(value >> precision) - precision) + 1;
}

class HLL{

private:
	int precision;
	int buckets;
	std::unique_ptr<std::vector<uint8_t> > table;
public:
	HLL(int precision)
	{
		this->precision = precision;
		this->buckets = (1ULL << precision);
		this->table =  std::make_unique<std::vector<uint8_t> > (this->buckets, 0);
    }

	~HLL() {}

	void add(uint32_t hash)
	{
		uint64_t bucket = hash & ((1ULL << this->precision) - 1);
		uint zeros = count_zeros (hash, this->precision) - 32;

		if (zeros > this->table->at (bucket))
		{
			this->table->at (bucket) = zeros;
		}
	}

	double get_threshold(int p)
	{
		return THRESHOLD_DATA[p - 4];
	}

	std::vector<int> get_nearest_neighbors(double E, std::vector<double> estimate)
	{
		std::vector< std::pair<double,int> > distance_map;
		std::vector<int> nearest;

		int i = 0;
		for (auto v: estimate) {
			std::pair<double, int> p(pow(E - v, 2.0), i);
			distance_map.push_back(p);
			i++;
		}

		sort(begin(distance_map), end(distance_map));

		for(int k=0; k < 6; k++) {
			nearest.push_back(distance_map[k].second);
		}

		return nearest;
	}

	double estimate_bias(double E, int p)
	{
		std::vector<int> nearest = get_nearest_neighbors(E, rawEstimateData[p]);
		double estimate = 0.0;

		for (auto v: nearest) {
			estimate += biasData[p][v];
		}
		return estimate / nearest.size();
	}


	size_t simpleQuery()
	{
		double sum = 0;

		for (size_t i = 0; i < this->buckets; ++i)
		{
			uint8_t bucket = this->table->at (i);
			sum += pow(2.0, -float(bucket));
		}

		double m = this->buckets;
		double harmonic = m / sum;

		double alpha;

		switch (this->precision) {
			case 4: alpha = 0.673; break;
			case 5: alpha = 0.697; break;
			case 6: alpha = 0.709; break;
			default:
				alpha = (0.7213 / (1.0 + (1.079 / double(1ULL << this->precision))));
		}

		long zeros = std::count(this->table->begin (), this->table->end (), 0);
		double estimate = alpha * m * harmonic;

		if (zeros > 0)
		{
			double H = m * std::log(m / float (zeros));
			if (H <= get_threshold(this->precision)) return H;
		}

		if (estimate <= (5 * m))
		{
			return estimate - estimate_bias (estimate, this->precision);
		}

		return estimate;
    }
};

#endif // HLL_SKETCH_H_INCLUDED
