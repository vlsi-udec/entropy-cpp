#include <cstdint>
#include <iostream>
#include <vector>
#include <array>
#include <functional>

#ifndef PACKET_H
#define PACKET_H

class packet {
public:
	int srcPort;
	int dstPort;
	int protocol;
	uint32_t srcIp;
	uint32_t dstIp;
	std::function<uint32_t (uint32_t *, uint32_t, uint32_t)> h_func;
	std::array<uint32_t, 5> temp;
	
	packet (int srcPort, int dstPort, uint32_t srcIp, uint32_t dstIp, int protocol,
	        std::function<uint32_t (uint32_t *, uint32_t, uint32_t)> h_func)
	{
		this->srcPort = srcPort;
		this->dstPort = dstPort;
		this->srcIp = srcIp;
		this->dstIp = dstIp;
		this->protocol = protocol;
		this->h_func = h_func;

		temp[0] = this->srcPort;
		temp[1] = this->dstPort;
		temp[2] = this->srcIp;
		temp[3] = this->dstIp;
		temp[4] = this->protocol;
	}

	~packet() {}

	uint64_t hash (uint32_t seed)
	{
		if (this->h_func) return this->h_func(temp.data(), seed, temp.size () * sizeof (uint32_t));
		return 0;
	}

	bool operator==(const packet& p) const
	{
		if (this->srcPort != p.srcPort) return false;
		if (this->dstPort != p.dstPort) return false;
		if (this->srcIp != p.srcIp) return false;
		if (this->dstIp != p.dstIp) return false;
		if (this->protocol != p.protocol) return false;
		return true;
	}
};

namespace std {
	template <>
	struct hash<packet>
	{
		std::size_t operator()(const packet& p) const
		{
			using std::size_t;
			using std::hash;
			using std::string;

			size_t spHash = std::hash<int>()(p.srcPort);
			size_t dpHash = std::hash<int>()(p.dstPort) << 1;
			size_t siHash = std::hash<int>()(p.srcIp) << 2;
			size_t diHash = std::hash<int>()(p.dstIp) << 3;
			size_t pHash = std::hash<int>()(p.protocol) << 4;
			return spHash ^ dpHash ^ siHash ^ diHash ^ pHash;
		}
	};
}

#endif // PACKET_H

