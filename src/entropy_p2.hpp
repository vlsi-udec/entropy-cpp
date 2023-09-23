#include <iostream>
#include <cmath>
#include "packet.hpp"
#include "entropy_base.hpp"
#include "pq_array.hpp"

#include "cms-cu.hpp"
#include "cms.hpp"
#include "cs.hpp"

#include "hll_sketch.h"

#ifndef ENTROPY_P2_H
#define ENTROPY_P2_H

template<class GenericCountSketch, class GenericPQ>
class entropy_p2: public entropy_base {
public:
	entropy_p2 (int w, int d, int queue_elements, uint q_bits, uint32_t seed, uint32_t seed_hll, int hll_p): cu(w, d), pq (queue_elements, q_bits, 36, seed), hll_sketch (hll_p), M(0), seed_hll(seed_hll) {
		queue_top_k = queue_elements * (1 << q_bits);
	}
	~entropy_p2 () {}

	void add_packet (packet p)
	{
		M++;
		int32_t est = cu.InsertSketch (p);
		pq.add (p, est);
		hll_sketch.add (p.hash(seed_hll));
	}

	std::pair<double, size_t> compute ()
	{
		size_t N_HLL = hll_sketch.simpleQuery ();
		auto hh_pq = pq.get_data ();
		uint32_t L_pq = 0;

		double entropy = 0;

		for (size_t i = 0; i < hh_pq.size(); i++)
		{
			L_pq += (hh_pq[i]);
			entropy += compute_entropy (hh_pq[i], M);
		}

		entropy -= compute_entropy_nhh (M, L_pq, queue_top_k, N_HLL);

		return std::make_pair (entropy, N_HLL);
	}

private:

	double compute_entropy_nhh (uint64_t M, uint32_t L, uint32_t K, uint32_t N)
	{
		return ((float)(M-L)/M)*log2((float)(M-L)/(M*(N-K)));
	}

	GenericCountSketch cu;
	GenericPQ pq;
	HLL hll_sketch;
	uint32_t seed_hll;
	size_t M;
	size_t queue_top_k;
};

#endif // ENTROPY_P2_H

