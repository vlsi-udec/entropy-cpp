#include <cstdint>
#include <iostream>
#include <vector>
#include <unordered_set>
#include <algorithm>
#include "pq_array.hpp"
#include "packet.hpp"

using namespace std;

pq_array::pq_array (uint elements_per_queue, uint array_bits, uint c_bits, uint seed)
{
	this->seed = seed;
	this->c_bits = c_bits;
	this->array_bits = array_bits;
	this->queue_len = elements_per_queue;

	for (uint i = 0; i < (1U << array_bits); ++i)
	{
		this->array.push_back (std::make_shared <std::vector<pq_element>> (this->queue_len+1, 0));
	}
}

pq_array::~pq_array ()
{
}

void pq_array::add (packet id, pq_count_t count)
{
	if (count < 0) return;
	// Hash
	uint32_t hash = id.hash(this->seed);
	uint32_t array_idx = hash & ((1 << this->array_bits) - 1);
	uint32_t int_id = (hash >> this->array_bits) & ((1UL << (72-c_bits))-1);

	// Overwrite
	bool overwrite = false;
	for (size_t i = 0; i < this->queue_len; ++i)
	{
		if (this->array[array_idx]->at(i).id == int_id)
		{
			this->array[array_idx]->at(i).count = count;
			overwrite = true;
		}
	}

	if (!overwrite)
		this->array[array_idx]->at(this->queue_len) = pq_element (id, int_id, count);

	sort (this->array[array_idx]->begin(), this->array[array_idx]->end(),
			std::greater<pq_element>());
	
	// if (this->array[array_idx].size () > queue_len)
		// this->array[array_idx].pop_back ();
}

std::vector <pq_count_t> pq_array::get_data ()
{
	std::vector <pq_count_t> data;

	for (auto & v: array)
	{
		for (size_t i = 0; i < this->queue_len; ++i)
		{
			data.push_back (v->at(i).count & ((1UL<<c_bits)-1) );
		}
	}

	if (verbose >= 2)
	{
		for (size_t i = 0; i < data.size(); ++i)
		{
			std::cout << i << " " << data[i] << "\n";
		}
	}

	return data;
}

std::unordered_set <packet> pq_array::get_id ()
{
	std::unordered_set <packet> data;

	for (auto & v: array)
	{
		// for (auto & e: v) data.insert (e.c_id);
	}

	return data;
}

void ppq::add (packet id, pq_count_t count)
{
	array[id] = count;
}
std::vector <pq_count_t> ppq::get_data ()
{
	std::vector <pq_count_t> data;

	for (auto & v: array)
	{
		data.push_back (v.second);
	}

	sort (data.begin(), data.end(), std::greater<pq_count_t>());
	data.resize(queue_len);

	if (verbose >= 2)
	{
		for (size_t i = 0; i < data.size(); ++i)
		{
			std::cout << i << " " << data[i] << "\n";
		}
	}

	return data;
}
