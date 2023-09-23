#include "entropy_p3.hpp"
#include <cmath>
#include <cstdint>
#include <cstdio>

extern int verbose;
#define printvar128(var) {printf ("%s: ", #var); std::cout << int64_t(var) << "\n";}
#define printvar(var) {printf ("%s: ", #var); std::cout << (var) << "\n";}

	template<size_t table_bits, size_t frac_bits>
	__int128 log2_hw(uint32_t data)
	{
		// return std::log2(data) * (1<<frac_bits);

		uint64_t d = data;
		int log2i = int(std::log2(data));
		if (log2i > table_bits)
			d = (d >> (log2i - table_bits)) & ((1 << table_bits) - 1);
		else if (log2i < table_bits) 
			d = (d << (table_bits - log2i)) & ((1 << table_bits) - 1);

		uint64_t base = (1 << frac_bits) * std::log2f(1+d / float(1 << table_bits));
		uint64_t next = (1 << frac_bits) * std::log2f(1+(d+1) / float(1 << table_bits));
		uint64_t log2f = base;

		uint64_t m = data;
		if (log2i > table_bits + 5)
		{
			m = (m >> (log2i - (table_bits + 5))) & ((1 << 5) - 1);
			log2f = base + (m * (next-base) >> 5);
		}
		else if (log2i > table_bits)
		{
			m = (m << ((table_bits + 5) - log2i)) & ((1 << 5) - 1);
			log2f = base + (m * (next-base) >> 5);
		}

		return (log2i << frac_bits) + log2f;
	}
	
	template<size_t T>
	__int128 div_sub (__int128 value)
	{
		// ap_int <T> ret;
		return (__int128(1)<<T) / value;
		// return ret;
	}


	template<size_t frac_bits>
	uint64_t newton_raphson (uint iterations, uint64_t b)
	{
		// return (1UL<<frac_bits)/float(b);
		uint64_t x0, xi;

		size_t i;
		for (i = 32; i >= 0; --i) if ((b >> i) & 0x1) break;

		x0 = 1UL << (frac_bits - i);

		for (i = 0; i < iterations; ++i)
		{
			xi = ((2UL << (frac_bits))) - (((b * x0)));
			xi = ((x0 * xi) >> (frac_bits));
			x0 = xi;
		}

		return xi;
	}

	template<size_t frac_bits, int log2_frac, int exp_frac>
	__int128 pow_hw (uint32_t base, int32_t exp)
	{
		// return (1<<frac_bits) * std::pow (base, exp / float (1<<18));

		__int128 iexp = exp;
		if (exp < 0) iexp = -exp;

		__int128 log2b = log2_hw<15,log2_frac>(base);
		__int128 m_exp_log = log2b * iexp; // QX.36
		__int128 m_exp_int = (m_exp_log >> (log2_frac + exp_frac));

		if (exp < 0)
		{
			__int128 m_exp_frac = (1<<frac_bits) * std::pow (2, (((m_exp_int + 1) << 36) - m_exp_log) / float(1UL<<36));
			return m_exp_frac >> (m_exp_int + 1);
		}
		else
		{
			__int128 m_exp_frac = (1<<frac_bits) * std::pow (2, (m_exp_log & ((1UL<<36) - 1)) / float(1UL<<36));
			m_exp_frac = m_exp_frac << m_exp_int;
			return m_exp_frac;
		}
	}


template<class GenericCountSketch, class GenericPQ>
void entropy_p3<GenericCountSketch, GenericPQ>::add_packet (packet p)
{
	M++;
	int32_t est = cu.InsertSketch (p);
	pq.add (p, est);
	hll_sketch.add (p.hash(seed_hll));
}

template<class GenericCountSketch, class GenericPQ>
std::pair<double, size_t> entropy_p3<GenericCountSketch, GenericPQ>::compute ()
{
	size_t N_HLL = hll_sketch.simpleQuery ();
	auto hh_pq = pq.get_data ();
	int queue_step = 1;
	
	double entropy = 0;
	uint64_t L_pq = 0;

	for (size_t i = 0; i < hh_pq.size(); i++)
	{
		L_pq += (hh_pq[i]);
	}

	double h1 = this->compute_h1 (hh_pq, M);

	if (queues_to_sort * queue_elements < hh_pq.size())
	{
		queue_step = hh_pq.size() / (queues_to_sort * queue_elements);
		hh_pq.resize(queues_to_sort * queue_elements);
	}

	// Sort pq
	if (sort_all)
	{
		std::sort(hh_pq.begin(), hh_pq.end(), std::greater<pq_count_t>());
	}
	else
	{
		std::vector<pq_count_t> q[queue_elements];
		for (size_t i = 0; i < queues_to_sort; ++i)
		{
			for (int j = 0; j < queue_elements; ++j)
			{
				q[j].push_back (hh_pq[queue_elements*i + j]);
			}
		}
		hh_pq.clear();

		for (int j = 0; j < queue_elements; ++j)
		{
			std::sort(q[j].begin(), q[j].end(), std::greater<pq_count_t>());
			for (auto &c: q[j])
			{
				hh_pq.push_back(c);
			}
		}
	}

	double h2 = this->compute_h2 (hh_pq, queue_step, M, L_pq, queue_top_k, N_HLL);

	std::cout << "H1 " << h1 << "\n";
	std::cout << "H2 " << h2 << "\n";

	return std::make_pair(h1 + h2, N_HLL);
}

// template<class GenericCountSketch>
// double entropy_p3_sw<GenericCountSketch>::compute_entropy_mls (double M, double L, double K, double N,
// 		                                                  double alpha, double alpha_i, double offset)
// {
// 	double C_i = ((1 - ((L + offset) / M)) * (alpha + 1)) / (std::pow (N, alpha + 1) - std::pow (K + 1, alpha + 1));
// 	double pis2 = alpha * C_i;
// 	double const_i = std::log2 (C_i) * (1 - ((L + offset) / M));
//
// 	double cont_1 = -std::log (N) / ((alpha_i-1) * std::pow (N, alpha_i - 1)) - 1 / (std::pow (alpha_i - 1, 2) * std::pow (N, alpha_i - 1));
// 	double cont_2 = std::log (K + 1) / ((alpha_i - 1) * std::pow (K + 1, alpha_i - 1))
// 				   + 1 / (std::pow (alpha_i - 1, 2) * std::pow (K + 1, alpha_i - 1));
// 	double cont = std::log2(std::exp(1)) * (cont_1 + cont_2);
//
// 	return pis2 * cont + const_i;
// }
//

template<class GenericCountSketch, class GenericPQ>
double entropy_p3_sw<GenericCountSketch, GenericPQ>::compute_entropy_mls (double M, double L, double K, double N,
		                                                  double alpha, double alpha_i, double offset)
{
	double invM = 1/M;
	double alpha_plus_one = alpha + 1;
	double alpha_i_minus_one = alpha_i - 1;
	double logN = std::log(N);
	double logK1 = std::log(K+1);


	if (verbose > 1)
	{
		printvar (M);
		printvar (L);
		printvar (K);
		printvar (N);
		printvar (alpha);
		printvar (alpha_i);
		printvar (invM);
		printvar (alpha_plus_one);
		printvar (alpha_i_minus_one);
		printvar (logN);
		printvar (logK1);
	}

	double LdM = (1 - ((L + offset) * invM));
	double C_i_num = (LdM * alpha_plus_one);
	double C_i_den = std::pow(N, alpha_plus_one) - std::pow(K+1, alpha_plus_one);
	double C_i = C_i_num / C_i_den;

	double pis2 = alpha  * C_i;
	double const_i = std::log2 (C_i) * LdM;

	double NPowAlpha_i = std::pow (N, alpha_i_minus_one);
	double K1PowAlpha_i = std::pow (K+1, alpha_i_minus_one);

	double cont_1_one = -logN / (alpha_i_minus_one * NPowAlpha_i);
	double cont_1_two = 1/(alpha_i_minus_one * alpha_i_minus_one * NPowAlpha_i);
	double cont_1 = cont_1_one - cont_1_two;

	double cont_2_one = logK1 / (alpha_i_minus_one * K1PowAlpha_i);
	double cont_2_two = 1/(alpha_i_minus_one * alpha_i_minus_one * K1PowAlpha_i);
	double cont_2 = cont_2_one + cont_2_two;

	double cont = std::log2(std::exp(1)) * (cont_1 + cont_2);

	if (verbose > 0)
	{
		std::cout << "LdM " << LdM << "\n";
		std::cout << "C_i_num " << C_i_num << "\n";
		std::cout << "C_i_den " << C_i_den << "\n";
		std::cout << "C_i " << C_i << "\n";

		std::cout << "pis2 " << pis2 << "\n";
		std::cout << "const_i " << const_i << "\n";

		std::cout << "NPowAlpha_i " << NPowAlpha_i << "\n";
		std::cout << "K1PowAlpha_i " << K1PowAlpha_i << "\n";

		std::cout << "cont_1_one_den " << alpha_i_minus_one * NPowAlpha_i << "\n";
		std::cout << "cont_1_one " << cont_1_one << "\n";
		std::cout << "cont_1_two " << cont_1_two << "\n";
		std::cout << "cont_1 " << cont_1 << "\n";

		std::cout << "cont_2_one " << cont_2_one << "\n";
		std::cout << "cont_2_two " << cont_2_two << "\n";
		std::cout << "cont_2 " << cont_2 << "\n";

		std::cout << "cont " << cont << "\n";
	}

	return pis2 * cont + const_i;
}


template<class GenericCountSketch, class GenericPQ>
double entropy_p3_sw<GenericCountSketch, GenericPQ>::compute_h1 (std::vector <pq_count_t> & hh_pq, size_t M)
{
	double entropy = 0;
	for (size_t i = 0; i < hh_pq.size(); i++)
	{
		entropy += entropy_base::compute_entropy (hh_pq[i], M);
	}
	return entropy;
}

template<class GenericCountSketch, class GenericPQ>
double entropy_p3_sw<GenericCountSketch, GenericPQ>::compute_h2 (std::vector <pq_count_t> & hh_pq, int queue_step,
		                                                      size_t M, size_t L, size_t K, size_t N)
{
	double entropy = 0;
	double sum_logi = 0;
	double sum_logmi = 0;
	double sum_logi2 = 0;
	double sum_logi_logmi = 0;

	for (size_t i = 0; i < hh_pq.size(); i++)
	{
		double logi = std::log2(queue_step*i+1);
		double logmi;

		if (hh_pq[i] == 0) logmi = 0;
		else logmi = std::log2(hh_pq[i]);

		sum_logi += logi;
		sum_logmi += logmi;
		sum_logi2 += logi * logi;
		sum_logi_logmi += logi * logmi;
	}

	double alpha_num = (hh_pq.size() * sum_logi_logmi - sum_logi * sum_logmi);
	double alpha_den = (hh_pq.size() * sum_logi2 - sum_logi * sum_logi);
	double alpha = alpha_num / alpha_den;
	double alpha_i = -alpha;
	double K2_f = (0 - std::log2(hh_pq[hh_pq.size() - 1])) / alpha + std::log2 (hh_pq.size());
	double K2 = std::round(std::pow(2, K2_f));

	if (verbose > 0)
	{
		printvar (alpha_num);
		printvar (alpha_den);
		printvar (alpha);
	}

	std::cout << "M " << M << "\n" ;
	std::cout << "N " << N << "\n" ;
	std::cout << "K2 " << K2 << "\n" ;

	if (K2 >= N)
	{
		// Two terms model
		entropy -= compute_entropy_mls (M, L, K, N, alpha, alpha_i);
	}
	else
	{
		// Three terms model
		entropy -= compute_entropy_mls (M, L, K, K2, alpha, alpha_i, N - K2);
		entropy += ((N - K2) / float (M)) * std::log2 (M);
	}
	return entropy;
}


template<class GenericCountSketch, class GenericPQ>
__int128 entropy_p3_hw<GenericCountSketch, GenericPQ>::compute_entropy_mls (__int128 M, __int128 L, __int128 K, __int128 N,
			                   __int128 alpha, __int128 alpha_i, __int128 offset)
{
	const uint64_t log2exp = (1<<log2_frac)*std::log2(std::exp(1));
	const uint64_t log2expinv = (1<<log2_frac)/std::log2(std::exp(1));

	__int128 invM = div_sub<div_frac> (M);
	__int128 alpha_plus_one = (alpha + (1UL << log2_frac));
	__int128 alpha_i_minus_one = (alpha_i - (1UL << log2_frac));
	__int128 logN = (log2_hw<log2_lut, log2_frac> (N) * log2expinv) >> log2_frac;// QX.18
	__int128 logK1 = (log2_hw<log2_lut, log2_frac> (K+1) * log2expinv) >> log2_frac;// QX.18

	__int128 LdM = (__int128(1)<<div_frac) - ((L + offset) * invM);

	__int128 C_i_num = (LdM * alpha_plus_one) >> div_frac; //QX.10
	__int128 C_i_den = pow_hw<pow_frac, log2_frac, log2_frac>(N, alpha_plus_one) - pow_hw<pow_frac, log2_frac, log2_frac>(K+1, alpha_plus_one);

	__int128 C_i = (C_i_num * div_sub<div_frac>(C_i_den)) >> (log2_frac); // QX.30
	__int128 pis2 = (alpha * C_i) >> log2_frac; // QX.30
	__int128 const_i = ((log2_hw<log2_lut, log2_frac> (C_i) - ((div_frac-pow_frac)<<log2_frac)) * LdM) >> log2_frac;

	__int128 NPowAlpha_i = pow_hw<pow_frac, log2_frac, log2_frac> (N, alpha_i_minus_one);
	__int128 K1PowAlpha_i = pow_hw<pow_frac, log2_frac, log2_frac> (K+1, alpha_i_minus_one);

	__int128 cont_1_one = (-logN * (div_sub<div_frac>((alpha_i_minus_one * NPowAlpha_i) >> log2_frac))) >> log2_frac;
	__int128 cont_1_two = div_sub<div_frac> ((alpha_i_minus_one * alpha_i_minus_one * NPowAlpha_i)>>(2*log2_frac));
	__int128 cont_1 = cont_1_one - cont_1_two;

	__int128 cont_2_one = (logK1 * (div_sub<div_frac>((alpha_i_minus_one * K1PowAlpha_i) >> log2_frac))) >> log2_frac;
	__int128 cont_2_two = div_sub<div_frac> ((alpha_i_minus_one * alpha_i_minus_one * K1PowAlpha_i)>>(2*log2_frac));
	__int128 cont_2 = cont_2_one + cont_2_two;

	__int128 cont = __int128(log2exp) * (cont_1 + cont_2); // QX.30
	__int128 ret = ((pis2*cont) >> (div_frac-pow_frac + log2_frac+div_frac - pow_frac - 2*log2_frac)) + (const_i >> (div_frac-2*log2_frac));

	if (verbose > 0)
	{
		std::cout << "logN " << float(logN)/float(1UL << log2_frac) << "\n";
		std::cout << "logK1 " << float(logK1)/float(1UL << log2_frac) << "\n";

		std::cout << "LdM " << float(LdM)/float(__int128(1) << div_frac) << "\n";
		std::cout << "C_i_num " << float(C_i_num)/float(1UL<< log2_frac) << "\n";
		std::cout << "C_i_den " << float(C_i_den)/float(1UL << (pow_frac)) << "\n";
		std::cout << "C_i " << float(C_i)/float(1UL<<(div_frac-pow_frac)) << "\n";

		std::cout << "pis2 " << float(pis2)/float(1UL<<(div_frac-pow_frac)) << "\n";
		std::cout << "const_i " << float(const_i)/float(__int128(1)<<div_frac) << "\n";
		//
		std::cout << "NPowAlpha_i " << float(NPowAlpha_i)/float(1UL<<pow_frac) << "\n";
		std::cout << "K1PowAlpha_i " << float(K1PowAlpha_i)/float(1UL<<pow_frac) << "\n";
		//
		std::cout << "cont_1_one_den " << float(alpha_i_minus_one * NPowAlpha_i)/float(1UL<<(pow_frac+log2_frac)) << "\n";
		std::cout << "cont_1_one " << float(cont_1_one)/float(1UL<<(div_frac - pow_frac)) << "\n";
		std::cout << "cont_1_two " << float(cont_1_two)/float(1UL<<(div_frac - pow_frac)) << "\n";
		std::cout << "cont_1 " << float(cont_1)/float(1UL<<(div_frac - pow_frac)) << "\n";

		std::cout << "cont_2_one " << float(cont_2_one)/float(1UL<<(div_frac - pow_frac)) << "\n";
		std::cout << "cont_2_two " << float(cont_2_two)/float(1UL<<(div_frac - pow_frac)) << "\n";
		std::cout << "cont_2 " << float(cont_2)/float(1UL<<(div_frac - pow_frac)) << "\n";

		std::cout << "cont " << float(cont)/float(1UL<<(log2_frac+div_frac - pow_frac)) << "\n";
		std::cout << "cont " << float(((pis2*cont) >> (div_frac-pow_frac +div_frac - pow_frac-log2_frac)))/float(1UL<<(2*log2_frac)) << "\n";
		
	}
	return ret;
}


// ap_int<60> entropy_p3_hw<GenericCountSketch>::compute_entropy_mls (uint64_t M, uint64_t L, uint64_t K,
//                                           __int128 N, int32_t alpha, int32_t alpha_i, uint64_t offset)
// {
	// const uint64_t log2exp = (1<<18)*std::log2(std::exp(1));
	// const uint64_t log2expinv = (1<<18)/std::log2(std::exp(1));

	// ap_int<32> invM = div_sub<div_frac> (M, 0);
	// ap_int<32> alpha_plus_one = (alpha + (1UL << 18));
	// ap_int<32> alpha_i_minus_one = (alpha_i - (1UL << 18));
	// ap_int<32> logN = (log2_hw<15, 18> (N) * log2expinv) >> 18;// QX.18
	// ap_int<32> logK1 = (log2_hw<15, 18> (K+1) * log2expinv) >> 18;// QX.18
	//
	// ap_int<32> LdM = ap_int<35>(1UL<<32) - (ap_int<35>(L + offset) * invM);
	//
	// ap_int<60> C_i_num = (LdM * alpha_plus_one) >> 30; //QX.10
	// ap_int<32> C_i_den = pow_hw<19>(N, alpha_plus_one.value) - pow_hw<19>(K+1, alpha_plus_one.value);
	// ap_int<50> C_i_den_inv = div_sub<div_frac>(C_i_den);
	//
	// __int128 C_i_p0 = (C_i_num.value & 0xFFFFFFFF) * (C_i_den_inv.value & 0xFFFFFFFF);
	// __int128 C_i_p1 = (C_i_num.value & 0xFFFFFFFF) * (C_i_den_inv.value >> 32);
	// __int128 C_i_p2 = (C_i_num.value >> 32) * (C_i_den_inv.value & 0xFFFFFFFF);
	// __int128 C_i_p3 = (C_i_num.value >> 32) * (C_i_den_inv.value >> 32);
	//
	// ap_int<60> C_i = (C_i_p0 >> 21) + (C_i_p1 << (32-21)) + (C_i_p2 << (32-21)) + (C_i_p3 << (64-21));
	// // ap_int<48> C_i = (C_i_num * div_sub<45>(C_i_den)) >> 16; // QX.30
	//
	// ap_int<48> pis2 = (ap_int<40> (alpha) * C_i) >> 18; // QX.30
	// 
	// ap_int<35> const_i = (ap_int<35>(((int32_t) log2_hw<15, 18> (C_i.value)) - (30<<18)) * LdM) >> 18;
	//
	// ap_int<40> NPowAlpha_i = pow_hw<30> (N, alpha_i_minus_one.value); // QX.25
	// ap_int<40> K1PowAlpha_i = pow_hw<19> (K+1, alpha_i_minus_one.value); // QX.19
	//
	// ap_int<50> cont_1_one_den = -div_sub<div_frac>((alpha_i_minus_one * NPowAlpha_i)>>22);
	// int64_t cont_1_one_p0 = (logN.value & 0xFFFFFFFF) * (cont_1_one_den.value & 0xFFFFFFFF);
	// int64_t cont_1_one_p1 = (logN.value & 0xFFFFFFFF) * (cont_1_one_den.value >> 32);
	// int64_t cont_1_one_p2 = (logN.value >> 32) * (cont_1_one_den.value & 0xFFFFFFFF);
	// int64_t cont_1_one_p3 = (logN.value >> 32) * (cont_1_one_den.value >> 32);
	//
	// ap_int<70> cont_1_one = (cont_1_one_p0 >> 18) + (cont_1_one_p1 << (32-18)) + (cont_1_one_p2 << (32-18)) + (cont_1_one_p3 << (64-18));
	//
	// ap_int<50> cont_1_two = div_sub<div_frac> ((((alpha_i_minus_one * alpha_i_minus_one)) * NPowAlpha_i)>>28);
	// ap_int<50> cont_1 = cont_1_one - cont_1_two;
	//
	//
	// ap_int<50> cont_2_one_den = div_sub<div_frac> (alpha_i_minus_one * K1PowAlpha_i >> 17);
	// int64_t cont_2_one_p0 = (logK1.value & 0xFFFFFFFF) * (cont_2_one_den.value & 0xFFFFFFFF);
	// int64_t cont_2_one_p1 = (logK1.value & 0xFFFFFFFF) * (cont_2_one_den.value >> 32);
	// int64_t cont_2_one_p2 = (logK1.value >> 32) * (cont_2_one_den.value & 0xFFFFFFFF);
	// int64_t cont_2_one_p3 = (logK1.value >> 32) * (cont_2_one_den.value >> 32);
	//
	// ap_int<50> cont_2_one = (cont_2_one_p0 >> 18) + (cont_2_one_p1 << (32-18)) + (cont_2_one_p2 << (32-18)) + (cont_2_one_p3 << (64-18));
	//
	// ap_int<40> cont_2_two = div_sub<div_frac> (alpha_i_minus_one * alpha_i_minus_one * K1PowAlpha_i >> 23);
	// ap_int<50> cont_2 = cont_2_one + cont_2_two;
	//
	// ap_int<64> cont = ap_int<32>(log2exp) * ((cont_1 + cont_2) >> 20); // QX.30
	// ap_int<60> ret = ((pis2*cont) >> 14) + (const_i << 4);
	//
	//
	// if (verbose > 0)
	// {
	// 	printvar(invM.to_fixed(32));
	// 	printvar(alpha_plus_one.to_fixed(18));
	// 	printvar128((LdM).value);
	// 	printvar128((alpha_plus_one).value);
	// 	printvar128((LdM * alpha_plus_one).value);
	//
	// 	std::cout << "LdM " << LdM.to_fixed(32) << "\n";
	// 	std::cout << "C_i_num " << C_i_num.to_fixed(20) << "\n";
	// 	std::cout << "C_i_den " << C_i_den.to_fixed(19) << "\n";
	// 	std::cout << "C_i " << C_i.to_fixed(30) << "\n";
	//
	// 	std::cout << "pis2 " << pis2.to_fixed(30) << "\n";
	// 	std::cout << "const_i " << const_i.to_fixed(32) << "\n";
	//
	// 	std::cout << "NPowAlpha_i " << NPowAlpha_i.to_fixed(30) << "\n";
	// 	std::cout << "K1PowAlpha_i " << K1PowAlpha_i.to_fixed(19) << "\n";
	//
	// 	std::cout << "cont_1_one_den " << cont_1_one_den.to_fixed(26) << "\n";
	// 	std::cout << "cont_1_one " << cont_1_one.to_fixed(22) << "\n";
	// 	std::cout << "cont_1_two " << cont_1_two.to_fixed(22) << "\n";
	// 	std::cout << "cont_1 " << cont_1.to_fixed(22) << "\n";
	//
	// 	std::cout << "cont_2_one " << cont_2_one.to_fixed(22) << "\n";
	// 	std::cout << "cont_2_two " << cont_2_two.to_fixed(22) << "\n";
	// 	std::cout << "cont_2 " << cont_2.to_fixed(22) << "\n";
	//
	// 	std::cout << "cont " << cont.to_fixed(30) << "\n";
	// }

	// return 0;
// }

template<class GenericCountSketch, class GenericPQ>
double entropy_p3_hw<GenericCountSketch, GenericPQ>::compute_h1 (std::vector <pq_count_t> & hh_pq, size_t M)
{
	uint64_t invM = newton_raphson<36> (15, M);
	uint64_t log2M = log2_hw<log2_lut, log2_frac> (M);

	double entropy = 0;
	int64_t h1 = 0;
	uint64_t L_pq = 0;
	for (size_t i = 0; i < hh_pq.size(); i++)
	{
		L_pq += (hh_pq[i]);
		h1 += hh_pq[i] * log2_hw<log2_lut, log2_frac> (hh_pq[i]);
	}
	h1 = ((L_pq * log2M - h1) * invM) >> (36-log2_frac); // Q28.36
	return h1 / float(1UL<<36);
}

template<class GenericCountSketch, class GenericPQ>
double entropy_p3_hw<GenericCountSketch, GenericPQ>::compute_h2 (std::vector <pq_count_t> & hh_pq, int queue_step,
		                                                      size_t M, size_t L, size_t K, size_t N)
{
	uint64_t invM = newton_raphson<36> (15, M);
	uint64_t log2M = log2_hw<log2_lut, log2_frac> (M);

	double sum_logi = 0;
	double sum_logi2 = 0;

	__int128 sum_logmi = 0;
	__int128 sum_logi_logmi = 0;

	for (size_t i = 0; i < hh_pq.size(); i++)
	{
		double logi_f = std::log2 (queue_step*i+1);
		__int128 logi = log2_hw<log2_lut, log2_frac> (queue_step*i+1);
		__int128 logmi = log2_hw<log2_lut, log2_frac> (hh_pq[i]);

		sum_logi += logi_f;
		sum_logi2 += logi_f * logi_f;
		sum_logmi += logmi;
		sum_logi_logmi += logi * logmi;
	}

	__int128 pq_size = hh_pq.size();
	__int128 alpha_num = (pq_size * sum_logi_logmi - __int128(sum_logi*(1<<log2_frac)) * sum_logmi) >> log2_frac; // QX.32
	double alpha_den = (hh_pq.size() * sum_logi2 - sum_logi * sum_logi);
	__int128 inv_alpha_den = (1UL<<div_frac) / alpha_den;

	__int128 alpha_p0 = (alpha_num >> 32) * (inv_alpha_den >> 32);
	__int128 alpha_p1 = (alpha_num >> 32) * (inv_alpha_den & 0xFFFFFFFF);
	__int128 alpha_p2 = (alpha_num & 0xFFFFFFFF) * (inv_alpha_den >> 32);
	__int128 alpha_p3 = (alpha_num & 0xFFFFFFFF) * (inv_alpha_den & 0xFFFFFFFF);

	__int128 alpha = (alpha_num * inv_alpha_den) >> div_frac;
	__int128 alpha_i = -alpha;

	__int128 K2_one = log2_hw<log2_lut, log2_frac> (hh_pq[hh_pq.size() - 1]);
	__int128 K2_two = (1<<log2_frac)*std::log2 (hh_pq.size());
	__int128 alpha_inv = div_sub<div_frac> (alpha_i); // Q1.31
	__int128 K2_i = ((alpha_inv * K2_one) >> (div_frac - log2_frac)) + K2_two;
	__int128 K2 = pow_hw<pow_frac, log2_frac, log2_frac> (2, K2_i) >> (pow_frac);
	
	double entropy;

	// if (verbose > 0)
	// {
	// 	printvar ((pq_size * sum_logi_logmi).to_fixed(32));
	// 	printvar ((ap_int<40>(sum_logi*(1<<14)) * sum_logmi).to_fixed(32));
	// 	printvar (alpha_num.to_fixed(32));
	// 	printvar(alpha_den);
	// 	printvar (alpha.to_fixed(18));
	// 	printvar (K2_i.to_fixed(18));
	// }
	std::cout << "N " << N << "\n" ;
	std::cout << "alpha_num " << float(alpha_num)/float(1<<log2_frac) << "\n" ;
	std::cout << "alpha_den " << float(alpha_den) << "\n" ;
	std::cout << "alpha_den_inv " << float(inv_alpha_den)/float(__int128(1)<<div_frac) << "\n" ;
	std::cout << "alpha " << float(alpha)/float(1<<log2_frac) << "\n" ;
	std::cout << "K2 " << float(K2) << "\n" ;

	if (K2 >= N)
	{
		// Two terms model
		__int128 h2 = compute_entropy_mls (M, L, K, N, alpha, alpha_i);
		entropy = -h2 / float(1UL<<36);
	}
	else
	{
		// Three terms model
		__int128 h2 = compute_entropy_mls (M, L, K, K2, alpha, alpha_i, N - K2);
		__int128 h3 = (((N - K2) *  invM) * log2M) >> (36-log2_frac);

		entropy = (h3 - h2) / float(1UL<<36);
	}
	return entropy;
}

template class entropy_p3_sw<CountSketch, pq_array>;
template class entropy_p3_sw<CountminSketch, pq_array>;
template class entropy_p3_sw<CountCuSketch, pq_array>;
template class entropy_p3_sw<CountPerfect, pq_array>;

template class entropy_p3_sw<CountSketch, ppq>;
template class entropy_p3_sw<CountminSketch, ppq>;
template class entropy_p3_sw<CountCuSketch, ppq>;
template class entropy_p3_sw<CountPerfect, ppq>;

template class entropy_p3_hw<CountSketch, pq_array>;
template class entropy_p3_hw<CountminSketch, pq_array>;
template class entropy_p3_hw<CountPerfect, pq_array>;
template class entropy_p3_hw<CountCuSketch, pq_array>;

template class entropy_p3_hw<CountSketch, ppq>;
template class entropy_p3_hw<CountminSketch, ppq>;
template class entropy_p3_hw<CountPerfect, ppq>;
template class entropy_p3_hw<CountCuSketch, ppq>;
