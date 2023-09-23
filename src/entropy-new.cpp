#include <cstdint>
#include <cstdio>
#include <memory>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <vector>
#include <cmath>
#include <algorithm>
#include <queue>
#include <set>
#include <utility>
#include <string>
#include <iomanip>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <limits>
#include <cstddef>
#include <omp.h>

#include "cms-cu.hpp"
#include "cms.hpp"
#include "cs.hpp"
#include "cm_perfect.hpp"
#include "hll_sketch.h"
#include "murmurhash.hpp"
#include "xxhash64.h"
#include "pq_array.hpp"
#include "packet.hpp"
#include "ap.hpp"

#include "entropy_p1.hpp"
#include "entropy_p2.hpp"
#include "entropy_p3.hpp"

#include "pcapplusplus/PcapFileDevice.h"
#include "pcapplusplus/Packet.h"
#include "pcapplusplus/RawPacket.h"
#include "pcapplusplus/IPv4Layer.h"
#include "pcapplusplus/TcpLayer.h"
#include "pcapplusplus/UdpLayer.h"

#include "malloc_count.h"

int verbose = 0;

static uint32_t xxhash_wrapper(const uint32_t *key, uint32_t seed, uint32_t len)
{
	return XXHash64::hash (key, len, seed);
}

int main(int argc, char *argv[])
{
	enum MODE {REAL = 0, P1 = 1, P2 = 2, P3 = 3, P3_HW = 4};

	char c;
	MODE mode = REAL;

	string filetype = "pcap";
	string filename = "";
	string hash = "murmurhash";
	string cm_type = "cs";
	string pq_type = "pqa";

	uint p_bits = 13;
	uint w_bits = 16;
	uint d_bits = 11;
	uint q_elements = 8;
	uint q_bits = 11;
	int queues_to_sort = -1;
	bool complete_queue_sort = true;

	while ((c = getopt(argc, argv, "t:f:h:p:w:d:e:q:m:s:k:lv:x:")) != -1)
	{
		switch (c) {
			case 't':
				filetype = std::string (optarg);
			break;
			case 'f':
				filename = std::string (optarg);
			break;
			case 's':
				cm_type = std::string (optarg);
			break;
			case 'x':
				pq_type = std::string (optarg);
			break;
			case 'h':
				hash = std::string (optarg);
			break;
			case 'p':
				p_bits = std::stoi (optarg);
			break;
			case 'w':
				w_bits = std::stoi (optarg);
			break;
			case 'd':
				d_bits = std::stoi (optarg);
			break;
			case 'e':
				q_elements = std::stoi (optarg);
			break;
			case 'q':
				q_bits = std::stoi (optarg);
			break;
			case 'm':
				mode = MODE(std::stoi (optarg));
			break;
			case 'k':
				queues_to_sort = std::stoi (optarg);
			break;
			case 'l':
				complete_queue_sort = false;
			break;
			case 'v':
				verbose = std::stoi (optarg);
			break;
			default:
			break;
		}
	}

	std::function<uint32_t (uint32_t *, uint32_t, uint32_t)> h_func;
	if (hash == "xxhash")
		h_func = xxhash_wrapper;
	else
		h_func = murmurhash;

	srand (10);
	uint32_t seed = rand ();
	uint32_t seed_hll = rand ();

	std::unique_ptr<entropy_base> entropy_handler;
	auto base_time = std::chrono::steady_clock::now ();
	size_t base_mem = malloc_count_current();

	switch (mode)
	{
		case MODE::REAL: {
			entropy_handler = std::make_unique<entropy_real> ();
		} break;
		case MODE::P1: {
			if (cm_type == "cs")
			{
				if (pq_type == "pqa")
				{
					entropy_handler = std::make_unique<entropy_p1<CountSketch, pq_array>> (w_bits, d_bits, q_elements, q_bits, seed);
				}
				else if (pq_type == "ppq")
				{
					entropy_handler = std::make_unique<entropy_p1<CountSketch, ppq>> (w_bits, d_bits, q_elements, q_bits, seed);
				}
			}
			else if (cm_type == "cm")
			{
				if (pq_type == "pqa")
				{
					entropy_handler = std::make_unique<entropy_p1<CountminSketch, pq_array>> (w_bits, d_bits, q_elements, q_bits, seed);
				}
				else if (pq_type == "ppq")
				{
					entropy_handler = std::make_unique<entropy_p1<CountminSketch, ppq>> (w_bits, d_bits, q_elements, q_bits, seed);
				}
			}
			else if (cm_type == "cu")
			{
				if (pq_type == "pqa")
				{
					entropy_handler = std::make_unique<entropy_p1<CountCuSketch, pq_array>> (w_bits, d_bits, q_elements, q_bits, seed);
				}
				else if (pq_type == "ppq")
				{
					entropy_handler = std::make_unique<entropy_p1<CountCuSketch, ppq>> (w_bits, d_bits, q_elements, q_bits, seed);
				}
			}
			else if (cm_type == "cp")
			{
				if (pq_type == "pqa")
				{
					entropy_handler = std::make_unique<entropy_p1<CountPerfect, pq_array>> (w_bits, d_bits, q_elements, q_bits, seed);
				}
				else if (pq_type == "ppq")
				{
					entropy_handler = std::make_unique<entropy_p1<CountPerfect, ppq>> (w_bits, d_bits, q_elements, q_bits, seed);
				}
			}
			else
			{
				std::cerr << "No valid sketch\n";
			}
		} break;
		case MODE::P2: {
			if (cm_type == "cs")
			{
				if (pq_type == "pqa")
				{
					entropy_handler = std::make_unique<entropy_p2<CountSketch, pq_array>> (w_bits, d_bits, q_elements, q_bits, seed, seed_hll, p_bits);
				}
				else if (pq_type == "ppq")
				{
					entropy_handler = std::make_unique<entropy_p2<CountSketch, ppq>> (w_bits, d_bits, q_elements, q_bits, seed, seed_hll, p_bits);
				}
			}
			else if (cm_type == "cm")
			{
				if (pq_type == "pqa")
				{
					entropy_handler = std::make_unique<entropy_p2<CountminSketch, pq_array>> (w_bits, d_bits, q_elements, q_bits, seed, seed_hll, p_bits);
				}
				else if (pq_type == "ppq")
				{
					entropy_handler = std::make_unique<entropy_p2<CountminSketch, ppq>> (w_bits, d_bits, q_elements, q_bits, seed, seed_hll, p_bits);
				}
			}
			else if (cm_type == "cu")
			{
				if (pq_type == "pqa")
				{
					entropy_handler = std::make_unique<entropy_p2<CountCuSketch, pq_array>> (w_bits, d_bits, q_elements, q_bits, seed, seed_hll, p_bits);
				}
				else if (pq_type == "ppq")
				{
					entropy_handler = std::make_unique<entropy_p2<CountCuSketch, ppq>> (w_bits, d_bits, q_elements, q_bits, seed, seed_hll, p_bits);
				}
			}
			else if (cm_type == "cp")
			{
				if (pq_type == "pqa")
				{
					entropy_handler = std::make_unique<entropy_p2<CountPerfect, pq_array>> (w_bits, d_bits, q_elements, q_bits, seed, seed_hll, p_bits);
				}
				else if (pq_type == "ppq")
				{
					entropy_handler = std::make_unique<entropy_p2<CountPerfect, ppq>> (w_bits, d_bits, q_elements, q_bits, seed, seed_hll, p_bits);
				}
			}
			else
			{
				std::cerr << "No valid sketch\n";
			}
		} break;
		case MODE::P3: {
			if (cm_type == "cs")
			{
				if (pq_type == "pqa")
				{
					entropy_handler = std::make_unique<entropy_p3_sw<CountSketch, pq_array>> (w_bits, d_bits, q_elements, q_bits, seed, seed_hll, p_bits, complete_queue_sort, queues_to_sort);
				}
				else if (pq_type == "ppq")
				{
					entropy_handler = std::make_unique<entropy_p3_sw<CountSketch, ppq>> (w_bits, d_bits, q_elements, q_bits, seed, seed_hll, p_bits, complete_queue_sort, queues_to_sort);
				}
			}
			else if (cm_type == "cm")
			{
				if (pq_type == "pqa")
				{
					entropy_handler = std::make_unique<entropy_p3_sw<CountminSketch, pq_array>> (w_bits, d_bits, q_elements, q_bits, seed, seed_hll, p_bits, complete_queue_sort, queues_to_sort);
				}
				else if (pq_type == "ppq")
				{
					entropy_handler = std::make_unique<entropy_p3_sw<CountminSketch, ppq>> (w_bits, d_bits, q_elements, q_bits, seed, seed_hll, p_bits, complete_queue_sort, queues_to_sort);
				}
			}
			else if (cm_type == "cu")
			{
				if (pq_type == "pqa")
				{
					entropy_handler = std::make_unique<entropy_p3_sw<CountCuSketch, pq_array>> (w_bits, d_bits, q_elements, q_bits, seed, seed_hll, p_bits, complete_queue_sort, queues_to_sort);
				}
				else if (pq_type == "ppq")
				{
					entropy_handler = std::make_unique<entropy_p3_sw<CountCuSketch, ppq>> (w_bits, d_bits, q_elements, q_bits, seed, seed_hll, p_bits, complete_queue_sort, queues_to_sort);
				}
			}
			else if (cm_type == "cp")
			{
				if (pq_type == "pqa")
				{
					entropy_handler = std::make_unique<entropy_p3_sw<CountPerfect, pq_array>> (w_bits, d_bits, q_elements, q_bits, seed, seed_hll, p_bits, complete_queue_sort, queues_to_sort);
				}
				else if (pq_type == "ppq")
				{
					entropy_handler = std::make_unique<entropy_p3_sw<CountPerfect, ppq>> (w_bits, d_bits, q_elements, q_bits, seed, seed_hll, p_bits, complete_queue_sort, queues_to_sort);
				}
			}
			else
			{
				std::cerr << "No valid sketch\n";
			}
		} break;
		case MODE::P3_HW: {
			if (pq_type == "pqa")
			{
				entropy_handler = std::make_unique<entropy_p3_hw<CountSketch, pq_array>> (w_bits, d_bits, q_elements, q_bits, seed, seed_hll, p_bits, complete_queue_sort, queues_to_sort);
			}
			else if (pq_type == "ppq")
			{
				entropy_handler = std::make_unique<entropy_p3_hw<CountSketch, ppq>> (w_bits, d_bits, q_elements, q_bits, seed, seed_hll, p_bits, complete_queue_sort, queues_to_sort);
			}
		} break;
		default:
			std::cerr << "No valid mode\n";
		return 0;
	}

	size_t unknown = 0;
	// Medicion memoria - antes de leer paquetes
	auto bp_time = std::chrono::steady_clock::now ();
	size_t bp_mem = malloc_count_current();
	if (filetype == "pcap")
	{
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

				pcpp::IPv4Address srcIP = ip_layer->getSrcIPv4Address();
				pcpp::IPv4Address dstIP = ip_layer->getDstIPv4Address();

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
				else 
				{
					unknown++;
				}

				packet temp (srcPort, dstPort, srcIP.toInt (), dstIP.toInt (), int(header->protocol), h_func);
				// Insercion a estructura
				entropy_handler->add_packet (temp);
			}
			else
			{
				unknown++;
			}
		}
		std::cout << "unknown " << unknown << "\n";
	}
	else if (filetype == "custom")
	{
		ifstream infile (filename);
		if (!infile.is_open())
		{
			printf("Error opening custom format file\n");
			return -1;
		}

		uint32_t _src_ip, _dst_ip;
		uint16_t _src_port, _dst_port;

		while (!infile.eof ())
		{
			infile >> _src_ip >> _src_port;
			infile >> _dst_ip >> _dst_port;

			packet temp (_src_port, _dst_port, _src_ip, _dst_ip, 0, h_func);
			entropy_handler->add_packet (temp);
		}

		infile.close();
	}
	else
	{
		std::cerr << "Non recognized format\n";
		return -1;
	}

	// Medicion tiempo - post
	auto pp_time = std::chrono::steady_clock::now ();
	size_t pp_mem = malloc_count_current();

	std::pair<double, size_t> e;
	e = entropy_handler->compute ();

	auto pc_time = std::chrono::steady_clock::now ();
	size_t pc_mem = malloc_count_current();

	std::cout << "Entropy: " << e.first << "\n";
	std::cout << "N: " << e.second << "\n";
	std::cout << "Entropy norm: " << (e.first/std::log2(e.second)) << "\n";
	std::cout << "Base memory: " << base_mem << " bytes\n";
	std::cout << "Pre process memory: " << bp_mem - base_mem << " bytes\n";
	std::cout << "Post process memory: " << pp_mem - bp_mem << " bytes\n";
	std::cout << "Post compute memory: " << pc_mem - pp_mem << " bytes\n";

	size_t pb_time_diff = std::chrono::duration_cast<std::chrono::microseconds>(pp_time - bp_time).count();
	size_t pp_time_diff = std::chrono::duration_cast<std::chrono::microseconds>(pc_time - pp_time).count();

	std::cout << "Process time: " << pb_time_diff << " ms\n";
	std::cout << "Compute time: " << pp_time_diff << " ms\n";

	return 0;
}

