#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <cmath>
#include <algorithm>
#include <queue>
#include <unordered_set>
#include <utility>
#include <string>
#include <iomanip>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <limits>
#include <cstddef>
#include <omp.h>

#include "pcapplusplus/PcapFileDevice.h"
#include "pcapplusplus/Packet.h"
#include "pcapplusplus/RawPacket.h"
#include "pcapplusplus/IPv4Layer.h"
#include "pcapplusplus/TcpLayer.h"
#include "pcapplusplus/UdpLayer.h"

using namespace std;

class packet {
public:
	int srcPort;
	int dstPort;
	int protocol;
	uint32_t srcIp;
	uint32_t dstIp;
	
	packet (int srcPort, int dstPort, uint32_t srcIp, uint32_t dstIp, int protocol)
	{
		this->srcPort = srcPort;
		this->dstPort = dstPort;
		this->srcIp = srcIp;
		this->dstIp = dstIp;
		this->protocol = protocol;
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

	struct HashFunction
	{
		size_t operator()(const packet& p) const
		{
			size_t spHash = std::hash<int>()(p.srcPort);
			size_t dpHash = std::hash<int>()(p.dstPort) << 1;
			size_t siHash = std::hash<int>()(p.srcIp) << 2;
			size_t diHash = std::hash<int>()(p.dstIp) << 3;
			size_t pHash = std::hash<int>()(p.protocol) << 4;
			return spHash ^ dpHash ^ siHash ^ diHash ^ pHash;
		}
	};
};

int main(int argc, char *argv[])
{
	if (argc != 2)
	{
		cerr << "Missing argument\n";
		cerr << "entropy-card [pcap file]\n";
		return -1;
	}

	unordered_set<packet,packet::HashFunction> input_data;
	unordered_set<int> input_srcPort;
	unordered_set<int> input_dstPort;
	unordered_set<uint32_t> input_srcIp;
	unordered_set<uint32_t> input_dstIp;
	unordered_set<int> input_proto;

	string filename = argv[1];


	pcpp::PcapFileReaderDevice reader (filename.c_str ());
	if (!reader.open())
	{
		printf("Error opening the pcap file\n");
		return 1;
	}

	pcpp::RawPacket rawPacket;
	while (reader.getNextPacket (rawPacket))
	{
		pcpp::Packet parsedPacket(&rawPacket);

		if (parsedPacket.isPacketOfType(pcpp::IPv4))
		{
			int srcPort, dstPort;

			pcpp::IPv4Layer * ip_layer = (pcpp::IPv4Layer *) parsedPacket.getLayerOfType (pcpp::IPv4);
		
			pcpp::IPv4Address srcIP = ip_layer->getSrcIpAddress();
			pcpp::IPv4Address dstIP = ip_layer->getDstIpAddress();
			
			pcpp::iphdr * header = ip_layer->getIPv4Header ();
			
			if (header->protocol == pcpp::PACKETPP_IPPROTO_TCP)
			{
				pcpp::TcpLayer * layer = (pcpp::TcpLayer *) parsedPacket.getLayerOfType (pcpp::TCP);
				if (layer == nullptr) continue;
				pcpp::tcphdr * tcp_header = layer->getTcpHeader ();
				
				srcPort = tcp_header->portSrc;
				dstPort = tcp_header->portDst;
			}
			else if (header->protocol == pcpp::PACKETPP_IPPROTO_UDP)
			{
				pcpp::UdpLayer * layer = (pcpp::UdpLayer *) parsedPacket.getLayerOfType (pcpp::UDP);
				if (layer == nullptr) continue;
				pcpp::udphdr * udp_header = layer->getUdpHeader ();
				
				srcPort = udp_header->portSrc;
				dstPort = udp_header->portDst;
			}
			
			packet temp (srcPort, dstPort, srcIP.toInt (), dstIP.toInt (), int(header->protocol));
			
			input_data.insert(temp);
			input_srcPort.insert(srcPort);
			input_dstPort.insert(dstPort);
			input_srcIp.insert(srcIP.toInt ());
			input_dstIp.insert(dstIP.toInt ());
			input_proto.insert(int(header->protocol));
		}
	}

	cout << "+ Database : " << filename << endl;
	cout << "+ Packet: " << input_data.size () << endl;
	cout << "+ SrcIp: " << input_srcIp.size () << endl;
	cout << "+ dstIP: " << input_dstIp.size () << endl;
	cout << "+ srcPort: " << input_srcPort.size () << endl;
	cout << "+ dstPort: " << input_dstPort.size () << endl;
	cout << "+ Protocol: " << input_proto.size () << endl;

	return 0;
}
