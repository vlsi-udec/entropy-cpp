#include <iostream>
#include <cmath>
#include "packet.hpp"
#include "ap.hpp"

#ifndef ENTROPY_BASE_H
#define ENTROPY_BASE_H

extern int verbose;
class entropy_base {
public:
	entropy_base () {}
	~entropy_base () {}
	virtual void add_packet (packet p) {}
	virtual std::pair<double, size_t> compute () {return std::make_pair(0, 0);}

	// Operations
	virtual float log2 (float value) {return std::log2f(value);}
	virtual float div (float dividend, float divisor) {return dividend/divisor;}

	double compute_entropy (double freq, double total)
	{
		if (freq == 0) return 0;
		return (freq/total)*std::log2(total/freq);
	}

	double norm_entropy (double entropy, double hh)
	{
		return (entropy)/std::log2(hh);
	}
};

class entropy_real: public entropy_base {
public:
	entropy_real () {}
	~entropy_real () {}

	void add_packet (packet p)
	{
		if (map_hh.find(p) != map_hh.end())
			map_hh[p] ++;
		else
			map_hh[p] = 1;
	}

	std::pair<double, size_t> compute ()
	{
		size_t M = 0;

		for (auto it = map_hh.begin(); it != map_hh.end(); ++it)
		{
			M = M + it->second;
		}

		if (verbose > 1)
		{
			size_t i = 0;
			for (auto it = map_hh.begin(); it != map_hh.end(); ++it)
			{
				std::cout << i << " " <<  it->second << "\n";
			}
			i++;
		}

		std::cout << "M " << M << "\n" ;

		double entropy = 0;
		
		for (auto it = map_hh.begin(); it != map_hh.end(); ++it)
		{
			entropy += compute_entropy (it->second, M);
		}
		return std::make_pair (entropy, map_hh.size());
	}

private:
	std::unordered_map<packet, uint32_t> map_hh;
};

#endif // ENTROPY_BASE_H

