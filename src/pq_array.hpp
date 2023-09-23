#include <cstdint>
#include <iostream>
#include <sys/types.h>
#include <unordered_map>
#include <vector>
#include <memory>
#include <unordered_set>
#include <algorithm>
#include "packet.hpp"
#include "common.hpp"

#ifndef PQ_ARRAY
#define PQ_ARRAY

class ppq {
private:
	std::unordered_map<packet, pq_count_t> array;
	uint queue_len;
public:
	ppq (uint elements_per_queue, uint array_bits, uint c_bits, uint seed): queue_len(elements_per_queue*(1<<array_bits)) {}
	~ppq (){}
	void add (packet id, pq_count_t count);
	std::vector <pq_count_t> get_data ();
};

class pq_element {
public:
	// packet c_id;
	uint32_t id;
	pq_count_t count;
	pq_element (int a)
	{
		id = 0;
		count = 0;
	}
	pq_element (packet c_id, uint32_t id, pq_count_t count) : id (id), count (count) {};
	friend bool operator> (const pq_element& o1, const pq_element& o2)
	{
		return o1.count > o2.count;
	}
};

class pq_array
{
private:

	uint32_t seed;
	uint array_bits;
	uint c_bits;
	uint queue_len;
	std::vector<std::shared_ptr<std::vector<pq_element>>> array;
public:
	pq_array (uint elements_per_queue, uint array_bits, uint c_bits, uint seed);
	~pq_array ();
	void add (packet id, pq_count_t count);
	std::vector <pq_count_t> get_data ();
	std::unordered_set <packet> get_id ();
};

#endif // PQ_ARRAY
