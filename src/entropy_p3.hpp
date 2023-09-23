#include <iostream>
#include <cmath>
#include "packet.hpp"
#include "entropy_base.hpp"
#include "pq_array.hpp"

#include "cms-cu.hpp"
#include "cms.hpp"
#include "cs.hpp"
#include "cm_perfect.hpp"

#include "hll_sketch.h"
#include "common.hpp"

#ifndef ENTROPY_P3_H
#define ENTROPY_P3_H

template<class GenericCountSketch, class GenericPQ>
class entropy_p3: public entropy_base {
public:
	entropy_p3 (int w, int d, int queue_elements, uint q_bits, uint32_t seed,
			uint32_t seed_hll, int hll_p, bool sort_all, size_t queues_to_sort):
		cu(w, d), pq (queue_elements, q_bits, 40, seed), hll_sketch (hll_p), M(0),
		seed_hll(seed_hll), sort_all(sort_all), queues_to_sort(queues_to_sort)
	{
		queue_top_k = queue_elements * (1 << q_bits);
		this->queue_elements = queue_elements;
	}
	~entropy_p3 () {}

	void add_packet (packet p);
	std::pair<double, size_t> compute ();

private:
	virtual double compute_h1 (std::vector <pq_count_t> & hh_pq, size_t M) {return 0;}
	virtual double compute_h2 (std::vector <pq_count_t> & hh_pq, int queue_step, size_t M, size_t L, size_t K, size_t N) {return 0;}

	GenericCountSketch cu;
	GenericPQ pq;
	HLL hll_sketch;
	uint32_t seed_hll;
	size_t M;
	size_t queue_top_k;
	int queue_elements;
	bool sort_all;
	size_t queues_to_sort;
};


template<class GenericCountSketch, class GenericPQ>
class entropy_p3_sw: public entropy_p3<GenericCountSketch, GenericPQ> {
public:
	entropy_p3_sw (int w, int d, int queue_elements, uint q_bits, uint32_t seed,
			uint32_t seed_hll, int hll_p, bool sort_all, size_t queues_to_sort):
		entropy_p3<GenericCountSketch, GenericPQ> (w, d, queue_elements, q_bits, seed,
				                 seed_hll, hll_p, sort_all, queues_to_sort)
	{}
	~entropy_p3_sw () {}

	double compute_h1 (std::vector <pq_count_t> & hh_pq, size_t M);
	double compute_h2 (std::vector <pq_count_t> & hh_pq, int queue_step, size_t M,
			                                       size_t L, size_t K, size_t N);

	double compute_entropy_mls (double M, double L, double K, double N, double alpha,
			                                      double alpha_i, double offset = 0);
};


template<class GenericCountSketch, class GenericPQ>
class entropy_p3_hw: public entropy_p3<GenericCountSketch, GenericPQ> {
public:
	static const int log2_lut = 15;
	static const int log2_frac = 18;
	static const int div_frac = 60;
	static const int pow_frac = 30;
	static const int alpha_frac = 30;
	

	entropy_p3_hw (int w, int d, int queue_elements, uint q_bits, uint32_t seed,
			uint32_t seed_hll, int hll_p, bool sort_all, size_t queues_to_sort):
		entropy_p3<GenericCountSketch, GenericPQ> (w, d, queue_elements, q_bits, seed,
				                 seed_hll, hll_p, sort_all, queues_to_sort)
	{}
	~entropy_p3_hw () {}

	double compute_h1 (std::vector <pq_count_t> & hh_pq, size_t M);
	double compute_h2 (std::vector <pq_count_t> & hh_pq, int queue_step, size_t M,
			                                       size_t L, size_t K, size_t N);

	__int128 compute_entropy_mls (__int128 M, __int128 L, __int128 K, __int128 N,
			                   __int128 alpha, __int128 alpha_i, __int128 offset = 0);
};


#endif // ENTROPY_P1_H

