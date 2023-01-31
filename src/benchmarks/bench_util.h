#ifndef ENCODINGS_ROOT_SRC_BENCH_MB_ENC_H_
#define ENCODINGS_ROOT_SRC_BENCH_MB_ENC_H_

#include "../encodings.h"
#include "../stratified_storage.h"
#include <cstdlib>

namespace enc::bench {

template<class T, lib::Encoding encoding>
inline T predicate_helper(T predicate) {
  if constexpr (encoding == lib::NONE) {
    return predicate;
  }
  if constexpr (encoding == lib::DFE) {
    return lib::ToDFE(predicate);
  }
  if constexpr (encoding == lib::EDFE) {
    return lib::ToEDFE(predicate);
  }
}

struct BenchmarkReportBlock {
  uint64_t result_acc = 0;
  uint64_t timer_ns_acc = 0;
};

inline uint64_t rng(uint64_t x) {
  // https://en.wikipedia.org/wiki/Xorshift
  x ^= x << 13;
  x ^= x >> 7;
  x ^= x << 17;
  return x;
}

}// namespace enc::bench

#endif//ENCODINGS_ROOT_SRC_BENCH_MB_ENC_H_
