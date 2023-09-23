## Proposed streaming algorithm to estimate empirical entropy of network flows

C++ implementation of the proposed streaming algorithm to estimate the empirical entropy on large datasets with sublinear memory usage. The algorithm uses sketches to estimate the frequency and cardinality of network flows during an observation interval. Also, the algorithm only stores the frequencies of the most frequent flows and uses them to estimate the rest of the frequencies by assuming a power-law distribution.

## Running the script entropy-paper

Example:
```
./entropy-paper -f <pcap-file> -m3 -d 13 -w 15 -e 8 -q10 -scu -xpqa

```
## Arguments

-f : Name of the input pcap file
-m : Running mode has three options: zero is for the actual entropy, 1 for the Top-K paper, 2 for the Top-K paper + uniform distribution, and 3 for the current paper using the Top-K and Power Law
-v2: This parameter allows prints of the contents of the perfect priority queue (ppq) or priority queue array (pqa)
-d : dimensions of the count sketch, numbers of rows
-w : dimensions of the count sketch, numbers of buckets in each row
-e : number of queue elements
-q : number of queues
-s : sketch type (cm, cs, cu)
-x : priority queue array (pqa) or  perfect priority queue(ppq)

## References
1. Top-K paper
2. Top-K paper + uniform distribution
3. Top-K and Power Law distribution
