#ifndef ENCODINGS_ROOT_SRC_BENCH_MICROBENCHMARKS_H_
#define ENCODINGS_ROOT_SRC_BENCH_MICROBENCHMARKS_H_

#include "../wrappers/taxi_wrapper.h"
#include "bench_util.h"
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <vector>

namespace enc::bench {

#define DEFAULT_MICROBENCHMARK_FETCH_ITERATIONS 1'000'000

template<class signed_T, class unsigned_T>
class MicrobenchmarkData {
 public:
  const std::string _experiment;
  const size_t _repeats;
  const signed_T _predicate_natural;
  const unsigned_T _predicate_dfe;
  const signed_T _predicate_edfe;
  bool _use_dfe;
  bool _use_edfe;
  std::vector<size_t> _idx;
  lib::BitStratifiedColumn _col_bit_natural;
  lib::BitStratifiedColumn _col_bit_dfe;
  lib::BitStratifiedColumn _col_bit_edfe;
  lib::ByteStratifiedColumn _col_byte_natural;
  lib::ByteStratifiedColumn _col_byte_dfe;
  lib::ByteStratifiedColumn _col_byte_edfe;
  const size_t _dataset_size;

  inline MicrobenchmarkData() = delete;
  inline MicrobenchmarkData(MicrobenchmarkData const &) = delete;

  inline MicrobenchmarkData(
      std::string &experiment,
      size_t repeats,
      std::vector<signed_T> dataset,
      signed_T predicate,
      int32_t bits_natural,
      int32_t bits_dfe,
      int32_t bits_edfe,
      size_t idx_length = DEFAULT_MICROBENCHMARK_FETCH_ITERATIONS,
      bool use_dfe = true,
      bool use_edfe = true,
      size_t idx_seed = -1) : _experiment(experiment),
                              _repeats(repeats),
                              _predicate_natural(predicate),
                              _predicate_dfe(use_dfe ? lib::ToDFE<unsigned_T>(predicate) : 0),
                              _predicate_edfe(use_edfe ? lib::ToEDFE<signed_T>(predicate) : 0),
                              _col_bit_natural(lib::BitStratifiedColumn(dataset.size(), bits_natural)),
                              _col_bit_dfe(lib::BitStratifiedColumn(dataset.size(), bits_dfe)),
                              _col_bit_edfe(lib::BitStratifiedColumn(dataset.size(), bits_edfe)),
                              _col_byte_natural(lib::ByteStratifiedColumn(dataset.size(), bits_natural)),
                              _col_byte_dfe(lib::ByteStratifiedColumn(dataset.size(), bits_dfe)),
                              _col_byte_edfe(lib::ByteStratifiedColumn(dataset.size(), bits_edfe)),
                              _use_dfe(use_dfe),
                              _use_edfe(use_edfe),
                              _dataset_size(dataset.size()) {
    std::vector<unsigned_T> dataset_dfe(dataset.size());
    std::vector<signed_T> dataset_edfe(dataset.size());
#pragma omp parallel for schedule(static) default(none) shared(dataset, dataset_dfe, dataset_edfe, _use_dfe, _use_edfe)
    for (size_t i = 0; i < _dataset_size; i++) {
      if (_use_dfe) {
        try {
          dataset_dfe[i] = lib::ToDFE<unsigned_T>((unsigned_T) dataset[i]);
        } catch (...) {
          _use_dfe = false;
        }
      }
      if (_use_edfe) {
        try {
          dataset_edfe[i] = lib::ToEDFE<signed_T>(dataset[i]);
        } catch (...) {
          _use_edfe = false;
        }
      }
    }
    _idx.reserve(idx_length);
    for (size_t i = 0; i < idx_length; i++) {
      idx_seed = rng(idx_seed);
      _idx.push_back(idx_seed % dataset.size());
    }
    lib::SetColumn(_col_bit_natural, dataset);
    lib::SetColumn(_col_byte_natural, dataset);
    if (use_dfe) {
      lib::SetColumn(_col_bit_dfe, dataset_dfe);
      lib::SetColumn(_col_byte_dfe, dataset_dfe);
    }
    if (use_edfe) {
      lib::SetColumn(_col_bit_edfe, dataset_edfe);
      lib::SetColumn(_col_byte_edfe, dataset_edfe);
    }
  }
};

template<lib::Encoding encoding, lib::BVModifier bv_update = lib::bvSet>
inline BenchmarkReportBlock bench_scan_bytestratified_reverse(
    const lib::ByteStratifiedColumn &dataset,
    int32_t predicate) {
  BenchmarkReportBlock out;
  auto bitvector = lib::Bitvector(dataset);
  const auto timer_start = std::chrono::high_resolution_clock::now();
  if constexpr (encoding == lib::NONE) {
    lib::LTConstantReverse<int32_t, encoding, bv_update>(dataset, bitvector, predicate);
  }
  if constexpr (encoding == lib::DFE) {
    lib::LTConstantReverse<int32_t, encoding, bv_update>(dataset, bitvector, predicate);
  }
  if constexpr (encoding == lib::EDFE) {
    lib::LTConstantReverse<int32_t, encoding, bv_update>(dataset, bitvector, predicate);
  }
  const auto timer_stop = std::chrono::high_resolution_clock::now();
  out.timer_ns_acc = std::chrono::duration_cast<std::chrono::nanoseconds>(
                         timer_stop - timer_start)
                         .count();
  out.result_acc = bitvector.popcnt();
  return out;
}

template<class signed_T, class unsigned_T, lib::BVModifier bv_update = lib::bvSet>
inline void bench_scan_bytestratified_reverse_wrapper(const MicrobenchmarkData<signed_T, unsigned_T> &data) {
  std::string header;
  header = "," + std::to_string(data._dataset_size) + "," + data._experiment + ",bench-scan-byte_s-reverse,natural,";
  for (int i = 0; i < data._repeats; i++) {
    BenchmarkReportBlock result = bench_scan_bytestratified_reverse<lib::NONE, bv_update>(data._col_byte_natural, data._predicate_natural);
    std::cout
        << std::to_string(i) + header + std::to_string(result.timer_ns_acc) + "," + std::to_string(result.result_acc) << std::endl;
  }
  header = "," + std::to_string(data._dataset_size) + "," + data._experiment + ",bench-scan-byte_s-reverse,dfe,";
  if (data._use_dfe) {
    for (int i = 0; i < data._repeats; i++) {
      BenchmarkReportBlock result = bench_scan_bytestratified_reverse<lib::DFE, bv_update>(data._col_byte_dfe, (signed_T) data._predicate_dfe);
      std::cout
          << std::to_string(i) + header + std::to_string(result.timer_ns_acc) + "," + std::to_string(result.result_acc) << std::endl;
    }
  }
  header = "," + std::to_string(data._dataset_size) + "," + data._experiment + ",bench-scan-byte_s-reverse,edfe,";
  if (data._use_edfe) {
    for (int i = 0; i < data._repeats; i++) {
      BenchmarkReportBlock result = bench_scan_bytestratified_reverse<lib::EDFE, bv_update>(data._col_byte_edfe, data._predicate_edfe);
      std::cout
          << std::to_string(i) + header + std::to_string(result.timer_ns_acc) + "," + std::to_string(result.result_acc) << std::endl;
    }
  }
}

template<lib::Encoding encoding, lib::BVModifier bv_update = lib::bvSet>
inline BenchmarkReportBlock bench_scan_bytestratified_forward(const lib::ByteStratifiedColumn &dataset, int32_t predicate) {
  BenchmarkReportBlock out;
  auto bitvector = lib::Bitvector(dataset);
  const auto timer_start = std::chrono::high_resolution_clock::now();
  if constexpr (encoding == lib::NONE) {
    lib::LTConstantForward<int32_t, encoding, bv_update>(dataset, bitvector, predicate);
  }
  if constexpr (encoding == lib::DFE) {
    lib::LTConstantForward<int32_t, encoding, bv_update>(dataset, bitvector, predicate);
  }
  if constexpr (encoding == lib::EDFE) {
    lib::LTConstantForward<int32_t, encoding, bv_update>(dataset, bitvector, predicate);
  }
  const auto timer_stop = std::chrono::high_resolution_clock::now();
  out.timer_ns_acc = std::chrono::duration_cast<std::chrono::nanoseconds>(
                         timer_stop - timer_start)
                         .count();
  out.result_acc = bitvector.popcnt();
  return out;
}

template<class signed_T, class unsigned_T, lib::BVModifier bv_update = lib::bvSet>
inline void bench_scan_bytestratified_forward_wrapper(const MicrobenchmarkData<signed_T, unsigned_T> &data) {
  std::string header;
  header = "," + std::to_string(data._dataset_size) + "," + data._experiment + ",bench-scan-byte_s-forward,natural,";
  for (int i = 0; i < data._repeats; i++) {
    BenchmarkReportBlock result = bench_scan_bytestratified_forward<lib::NONE, bv_update>(data._col_byte_natural, data._predicate_natural);
    std::cout
        << std::to_string(i) + header + std::to_string(result.timer_ns_acc) + "," + std::to_string(result.result_acc) << std::endl;
  }
  header = "," + std::to_string(data._dataset_size) + "," + data._experiment + ",bench-scan-byte_s-forward,dfe,";
  if (data._use_dfe) {
    for (int i = 0; i < data._repeats; i++) {
      BenchmarkReportBlock result = bench_scan_bytestratified_forward<lib::DFE, bv_update>(data._col_byte_dfe, (signed_T) data._predicate_dfe);
      std::cout
          << std::to_string(i) + header + std::to_string(result.timer_ns_acc) + "," + std::to_string(result.result_acc) << std::endl;
    }
  }
  header = "," + std::to_string(data._dataset_size) + "," + data._experiment + ",bench-scan-byte_s-forward,edfe,";
  if (data._use_edfe) {
    for (int i = 0; i < data._repeats; i++) {
      BenchmarkReportBlock result = bench_scan_bytestratified_forward<lib::EDFE, bv_update>(data._col_byte_edfe, data._predicate_edfe);
      std::cout
          << std::to_string(i) + header + std::to_string(result.timer_ns_acc) + "," + std::to_string(result.result_acc) << std::endl;
    }
  }
}

template<lib::Encoding encoding>
inline BenchmarkReportBlock bench_fetch_bytestratified(const std::vector<size_t> &idx_list, const lib::ByteStratifiedColumn &dataset) {
  BenchmarkReportBlock out;
  for (const unsigned long long idx : idx_list) {
    uint64_t result;
    const auto timer_start = std::chrono::high_resolution_clock::now();
    if constexpr (encoding == lib::NONE) {
      result = lib::SimpleFetchFromColumn<int32_t>(dataset, idx);
    }
    if constexpr (encoding == lib::DFE) {
      result = lib::DfeFetchFromColumn<uint32_t>(dataset, idx);
    }
    if constexpr (encoding == lib::EDFE) {
      result = lib::EdfeFetchFromColumn<int32_t>(dataset, idx);
    }
    const auto timer_stop = std::chrono::high_resolution_clock::now();
    out.result_acc += result;
    out.timer_ns_acc += std::chrono::duration_cast<std::chrono::nanoseconds>(
                            timer_stop - timer_start)
                            .count();
  }
  return out;
}

template<class signed_T, class unsigned_T>
inline void bench_fetch_bytestratified_wrapper(const MicrobenchmarkData<signed_T, unsigned_T> &data) {
  std::string header;
  header = "," + std::to_string(data._dataset_size) + "," + data._experiment + ",bench-fetch-byte_s,natural,";
  for (int i = 0; i < data._repeats; i++) {
    BenchmarkReportBlock result = bench_fetch_bytestratified<lib::NONE>(data._idx, data._col_byte_natural);
    std::cout
        << std::to_string(i) + header + std::to_string(result.timer_ns_acc) + "," + std::to_string(result.result_acc) << std::endl;
  }
  header = "," + std::to_string(data._dataset_size) + "," + data._experiment + ",bench-fetch-byte_s,dfe,";
  if (data._use_dfe) {
    for (int i = 0; i < data._repeats; i++) {
      BenchmarkReportBlock result = bench_fetch_bytestratified<lib::DFE>(data._idx, data._col_byte_dfe);
      std::cout
          << std::to_string(i) + header + std::to_string(result.timer_ns_acc) + "," + std::to_string(result.result_acc) << std::endl;
    }
  }
  header = "," + std::to_string(data._dataset_size) + "," + data._experiment + ",bench-fetch-byte_s,edfe,";
  if (data._use_edfe) {
    for (int i = 0; i < data._repeats; i++) {
      BenchmarkReportBlock result = bench_fetch_bytestratified<lib::EDFE>(data._idx, data._col_byte_edfe);
      std::cout
          << std::to_string(i) + header + std::to_string(result.timer_ns_acc) + "," + std::to_string(result.result_acc) << std::endl;
    }
  }
}

template<lib::Encoding encoding, lib::BVModifier bv_update = lib::bvSet>
inline BenchmarkReportBlock bench_scan_bitstratified_forward(const lib::BitStratifiedColumn &dataset, int32_t predicate) {
  BenchmarkReportBlock out;
  auto bitvector = lib::Bitvector(dataset);
  const auto timer_start = std::chrono::high_resolution_clock::now();
  if constexpr (encoding == lib::NONE) {
    lib::LTConstantForward<int32_t, encoding, bv_update>(dataset, bitvector, predicate);
  }
  if constexpr (encoding == lib::DFE) {
    lib::LTConstantForward<int32_t, encoding, bv_update>(dataset, bitvector, predicate);
  }
  if constexpr (encoding == lib::EDFE) {
    lib::LTConstantForward<int32_t, encoding, bv_update>(dataset, bitvector, predicate);
  }
  const auto timer_stop = std::chrono::high_resolution_clock::now();
  out.timer_ns_acc = std::chrono::duration_cast<std::chrono::nanoseconds>(
                         timer_stop - timer_start)
                         .count();
  out.result_acc = bitvector.popcnt();
  return out;
}

template<class signed_T, class unsigned_T, lib::BVModifier bv_update = lib::bvSet>
inline void bench_scan_bitstratified_forward_wrapper(const MicrobenchmarkData<signed_T, unsigned_T> &data) {
  std::string header;
  header = "," + std::to_string(data._dataset_size) + "," + data._experiment + ",bench-scan-bit_s-forward,natural,";
  for (int i = 0; i < data._repeats; i++) {
    BenchmarkReportBlock result = bench_scan_bitstratified_forward<lib::NONE, bv_update>(data._col_bit_natural, data._predicate_natural);
    std::cout
        << std::to_string(i) + header + std::to_string(result.timer_ns_acc) + "," + std::to_string(result.result_acc) << std::endl;
  }
  header = "," + std::to_string(data._dataset_size) + "," + data._experiment + ",bench-scan-bit_s-forward,dfe,";
  if (data._use_dfe) {
    for (int i = 0; i < data._repeats; i++) {
      BenchmarkReportBlock result = bench_scan_bitstratified_forward<lib::DFE, bv_update>(data._col_bit_dfe, (signed_T) data._predicate_dfe);
      std::cout
          << std::to_string(i) + header + std::to_string(result.timer_ns_acc) + "," + std::to_string(result.result_acc) << std::endl;
    }
  }
  header = "," + std::to_string(data._dataset_size) + "," + data._experiment + ",bench-scan-bit_s-forward,edfe,";
  if (data._use_edfe) {
    for (int i = 0; i < data._repeats; i++) {
      BenchmarkReportBlock result = bench_scan_bitstratified_forward<lib::EDFE, bv_update>(data._col_bit_edfe, data._predicate_edfe);
      std::cout
          << std::to_string(i) + header + std::to_string(result.timer_ns_acc) + "," + std::to_string(result.result_acc) << std::endl;
    }
  }
}

template<lib::Encoding encoding, lib::BVModifier bv_update = lib::bvSet>
inline BenchmarkReportBlock bench_scan_bitstratified_reverse(const lib::BitStratifiedColumn &dataset, int32_t predicate) {
  BenchmarkReportBlock out;
  auto bitvector = lib::Bitvector(dataset);
  const auto timer_start = std::chrono::high_resolution_clock::now();
  if constexpr (encoding == lib::NONE) {
    lib::LTConstantReverse<int32_t, lib::NONE, bv_update>(dataset, bitvector, predicate);
  }
  if constexpr (encoding == lib::DFE) {
    lib::LTConstantReverse<int32_t, lib::DFE, bv_update>(dataset, bitvector, predicate);
  }
  if constexpr (encoding == lib::EDFE) {
    lib::LTConstantReverse<int32_t, lib::EDFE, bv_update>(dataset, bitvector, predicate);
  }
  const auto timer_stop = std::chrono::high_resolution_clock::now();
  out.timer_ns_acc = std::chrono::duration_cast<std::chrono::nanoseconds>(
                         timer_stop - timer_start)
                         .count();
  out.result_acc = bitvector.popcnt();
  return out;
}

template<class signed_T, class unsigned_T, lib::BVModifier bv_update = lib::bvSet>
inline void bench_scan_bitstratified_reverse_wrapper(const MicrobenchmarkData<signed_T, unsigned_T> &data) {
  std::string header;
  header = "," + std::to_string(data._dataset_size) + "," + data._experiment + ",bench-scan-bit_s-reverse,natural,";
  for (int i = 0; i < data._repeats; i++) {
    BenchmarkReportBlock result = bench_scan_bitstratified_reverse<lib::NONE, bv_update>(data._col_bit_natural, data._predicate_natural);
    std::cout
        << std::to_string(i) + header + std::to_string(result.timer_ns_acc) + "," + std::to_string(result.result_acc) << std::endl;
  }
  header = "," + std::to_string(data._dataset_size) + "," + data._experiment + ",bench-scan-bit_s-reverse,dfe,";
  if (data._use_dfe) {
    for (int i = 0; i < data._repeats; i++) {
      BenchmarkReportBlock result = bench_scan_bitstratified_reverse<lib::DFE, bv_update>(data._col_bit_dfe, (signed_T) data._predicate_dfe);
      std::cout
          << std::to_string(i) + header + std::to_string(result.timer_ns_acc) + "," + std::to_string(result.result_acc) << std::endl;
    }
  }
  header = "," + std::to_string(data._dataset_size) + "," + data._experiment + ",bench-scan-bit_s-reverse,edfe,";
  if (data._use_edfe) {
    for (int i = 0; i < data._repeats; i++) {
      BenchmarkReportBlock result = bench_scan_bitstratified_reverse<lib::EDFE, bv_update>(data._col_bit_edfe, data._predicate_edfe);
      std::cout
          << std::to_string(i) + header + std::to_string(result.timer_ns_acc) + "," + std::to_string(result.result_acc) << std::endl;
    }
  }
}

template<lib::Encoding encoding>
inline BenchmarkReportBlock bench_fetch_bitstratified(const std::vector<size_t> &idx_list, const lib::BitStratifiedColumn &dataset) {
  BenchmarkReportBlock out;
  for (const unsigned long long idx : idx_list) {
    uint64_t result;
    const auto timer_start = std::chrono::high_resolution_clock::now();
    if constexpr (encoding == lib::NONE) {
      result = lib::SimpleFetchFromColumn<int32_t>(dataset, idx);
    }
    if constexpr (encoding == lib::DFE) {
      result = lib::DfeFetchFromColumn<uint32_t>(dataset, idx);
    }
    if constexpr (encoding == lib::EDFE) {
      result = lib::EdfeFetchFromColumn<int32_t>(dataset, idx);
    }
    const auto timer_stop = std::chrono::high_resolution_clock::now();
    out.result_acc += result;
    out.timer_ns_acc += std::chrono::duration_cast<std::chrono::nanoseconds>(
                            timer_stop - timer_start)
                            .count();
  }
  return out;
}

template<class signed_T, class unsigned_T>
inline void bench_fetch_bitstratified_wrapper(const MicrobenchmarkData<signed_T, unsigned_T> &data) {
  std::string header;
  header = "," + std::to_string(data._dataset_size) + "," + data._experiment + ",bench-fetch-bit_s,natural,";
  for (int i = 0; i < data._repeats; i++) {
    BenchmarkReportBlock result = bench_fetch_bitstratified<lib::NONE>(data._idx, data._col_bit_natural);
    std::cout
        << std::to_string(i) + header + std::to_string(result.timer_ns_acc) + "," + std::to_string(result.result_acc) << std::endl;
  }
  header = "," + std::to_string(data._dataset_size) + "," + data._experiment + ",bench-fetch-bit_s,dfe,";
  if (data._use_dfe) {
    for (int i = 0; i < data._repeats; i++) {
      BenchmarkReportBlock result = bench_fetch_bitstratified<lib::DFE>(data._idx, data._col_bit_dfe);
      std::cout
          << std::to_string(i) + header + std::to_string(result.timer_ns_acc) + "," + std::to_string(result.result_acc) << std::endl;
    }
  }
  header = "," + std::to_string(data._dataset_size) + "," + data._experiment + ",bench-fetch-bit_s,edfe,";
  if (data._use_edfe) {
    for (int i = 0; i < data._repeats; i++) {
      BenchmarkReportBlock result = bench_fetch_bitstratified<lib::EDFE>(data._idx, data._col_bit_edfe);
      std::cout
          << std::to_string(i) + header + std::to_string(result.timer_ns_acc) + "," + std::to_string(result.result_acc) << std::endl;
    }
  }
}

template<class signed_T, class unsigned_T, lib::BVModifier bv_update = lib::bvSet>
inline void microbenchmark(const MicrobenchmarkData<signed_T, unsigned_T> &data) {
  bench_scan_bitstratified_forward_wrapper<signed_T, unsigned_T, bv_update>(data);
  bench_scan_bitstratified_reverse_wrapper<signed_T, unsigned_T, bv_update>(data);
  bench_scan_bytestratified_forward_wrapper<signed_T, unsigned_T, bv_update>(data);
  bench_scan_bytestratified_reverse_wrapper<signed_T, unsigned_T, bv_update>(data);
  bench_fetch_bitstratified_wrapper<signed_T, unsigned_T>(data);
  bench_fetch_bytestratified_wrapper<signed_T, unsigned_T>(data);
}

inline void microbenchmark_uniform(size_t size, size_t repeats) {
  constexpr uint32_t bench_mod_val = 10'000'000;
  constexpr int32_t predicate = 1'000'000;
  constexpr int32_t col_bits = 32;
  std::string experiment = "microbenchmark-uniform";
  uint64_t seed = -1;
  std::vector<int32_t> dataset;
  dataset.reserve(size);
  for (int i = 0; i < size; i++) {
    seed = rng(seed);
    int32_t this_val = ((uint32_t) seed) % bench_mod_val;
    dataset.push_back(this_val);
  }
  auto data = MicrobenchmarkData<int32_t, uint32_t>(
      experiment,
      repeats,
      dataset,
      predicate,
      col_bits,
      col_bits,
      col_bits);
  microbenchmark<int32_t, uint32_t>(data);
}

inline void microbenchmark_onezero(size_t size, size_t repeats) {
  constexpr int32_t predicate = 1;
  constexpr int32_t col_bits = 32;
  std::string experiment = "microbenchmark-onezero";
  uint64_t seed = -1;
  std::vector<int32_t> dataset;
  dataset.reserve(size);
  for (int i = 0; i < size; i++) {
    seed = rng(seed);
    int32_t this_val = ((uint32_t) seed) & 0x1;
    dataset.push_back(this_val);
  }
  auto data = MicrobenchmarkData<int32_t, uint32_t>(
      experiment,
      repeats,
      dataset,
      predicate,
      col_bits,
      col_bits,
      col_bits);
  microbenchmark<int32_t, uint32_t>(data);
}

inline void microbenchmark_taxi_trip_distance(const lib::TaxiTable &table, size_t repeats) {
  constexpr int32_t predicate = 100;
  constexpr int32_t col_bits = 32;
  std::string experiment = "microbenchmark-taxi-trip-distance";
  auto data = MicrobenchmarkData<int32_t, uint32_t>(
      experiment,
      repeats,
      table.trip_distance,
      predicate,
      col_bits,
      col_bits,
      col_bits);
  microbenchmark<int32_t, uint32_t, lib::bvSetNot>(data);
}

inline void microbenchmark_taxi_passenger_count(const lib::TaxiTable &table, size_t repeats) {
  constexpr int32_t predicate = 2;
  constexpr int32_t col_bits = 32;
  std::string experiment = "microbenchmark-taxi-passenger-count";
  auto data = MicrobenchmarkData<int32_t, uint32_t>(
      experiment,
      repeats,
      table.passenger_count,
      predicate,
      col_bits,
      col_bits,
      col_bits);
  microbenchmark<int32_t, uint32_t, lib::bvSetNot>(data);
}

inline void microbenchmark_taxi_tip_amount(const lib::TaxiTable &table, size_t repeats) {
  constexpr int32_t predicate = 1000;
  constexpr int32_t col_bits = 32;
  std::string experiment = "microbenchmark-taxi-tip-amount";
  auto data = MicrobenchmarkData<int32_t, uint32_t>(
      experiment,
      repeats,
      table.tip_amount,
      predicate,
      col_bits,
      col_bits,
      col_bits);
  microbenchmark<int32_t, uint32_t, lib::bvSetNot>(data);
}

inline void microbenchmarks(const lib::TaxiTable &table, size_t repeats) {
  microbenchmark_taxi_passenger_count(table, repeats);
  microbenchmark_taxi_tip_amount(table, repeats);
  microbenchmark_taxi_trip_distance(table, repeats);
  microbenchmark_onezero(table.n(), repeats);
  microbenchmark_uniform(table.n(), repeats);
}

}// namespace enc::bench

#endif//ENCODINGS_ROOT_SRC_BENCH_MICROBENCHMARKS_H_
