#include <iostream>
#include <errno.h>
#include <string>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <algorithm>
#include <unordered_set>
#include <bits/stdc++.h>

#ifndef PARAMETROS_H_INCLUDED
#define PARAMETROS_H_INCLUDED

inline double Alphas(int precision) {
  switch (precision) {
    case 4:
      return 0.673;
    case 5:
      return 0.697;
    case 6:
      return 0.709;
    default:
      return (0.7213 / (1.0 + (1.079 / static_cast<double>(1 << precision))));
  }
}

//get the index of M
inline int getIndexM(uint64_t hash, int precision) {
  return (hash & ((1 << precision)-1));
}

// count leading zeros
inline uint8_t countZeroes(uint64_t hash, int precision) {
  // const uint64_t mask = ~(((1 << precision) - (uint64_t)1) << (64 - precision));

  return (__builtin_clzll(hash >> precision) - static_cast<uint8_t>(precision));
}


#endif // PARAMETROS_H_INCLUDED
