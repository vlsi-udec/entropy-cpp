#include <iostream>
#include <cmath>
#include "packet.hpp"
#include "entropy_base.hpp"
#include "pq_array.hpp"

#include "cms-cu.hpp"
#include "cms.hpp"
#include "cs.hpp"

#ifndef ENTROPY_P1_H
#define ENTROPY_P1_H

template<class GenericCountSketch, class GenericPQ>
class entropy_p1: public entropy_base {
public:
	entropy_p1 (int w, int d, int queue_elements, uint q_bits, int seed): cu(w, d), pq (queue_elements, q_bits, 36, seed) {}
	~entropy_p1 () {}

	void add_packet (packet p)
	{
		int32_t est = cu.InsertSketch (p);
		pq.add (p, est);
	}

	std::pair<double, size_t> compute ()
	{
		auto hh_pq = pq.get_data ();
		uint32_t L_pq = 0;
		
		for (size_t i = 0; i < hh_pq.size(); i++)
		{
			L_pq += (hh_pq[i]);
		}

		double entropy = 0;
		for (size_t i = 0; i < hh_pq.size(); i++)
		{
			entropy += compute_entropy (hh_pq[i], L_pq);
		}

		return std::make_pair (entropy, hh_pq.size());
	}

private:
	GenericCountSketch cu;
	GenericPQ pq;
};

#endif // ENTROPY_P1_H

