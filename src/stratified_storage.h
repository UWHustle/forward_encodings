#ifndef ENCODINGS_ROOT_SRC_STRATIFIED_STORAGE_H_
#define ENCODINGS_ROOT_SRC_STRATIFIED_STORAGE_H_

#include <bitset>
#include <cstdlib>
#include <immintrin.h>
#include <stdexcept>
#include <vector>

#define BITVECTOR_PREFETCH_DISTANCE 1024
#define BIT_FETCH_ENABLE_PREFETCH true
#define BIT_FETCH_MINIMIZE_READS true
#define BYTE_FETCH_ENABLE_PREFETCH true
#define BYTE_FETCH_MINIMIZE_READS true
#define SCAN_ENABLE_PREFETCH true
#define SCAN_VERTICAL_PREFETCH_DISTANCE 1024

namespace enc::lib {

class ByteStratifiedColumn {
 public:
  const uint64_t val_bits;
  const uint64_t num_strata;
  const uint64_t right_padding;
  const size_t size;
  const size_t num_avx2_blocks;
  const size_t page_aligned_size;

  uint8_t **strata = nullptr;

  inline ByteStratifiedColumn() = delete;
  inline ByteStratifiedColumn(ByteStratifiedColumn const &) = delete;
  void operator=(ByteStratifiedColumn const &) = delete;

  inline ByteStratifiedColumn(
      size_t column_size,
      uint64_t value_bits) : val_bits((int) value_bits),
                             num_strata((int) (((value_bits % 8 != 0) ? 1 : 0) + (value_bits / 8))),
                             right_padding((num_strata * 8) - val_bits),
                             size(column_size),
                             num_avx2_blocks((column_size / 32) + ((column_size % 32) ? 1 : 0)),
                             page_aligned_size(4096 - (column_size % 4096) + column_size) {
    strata = static_cast<uint8_t **>(std::malloc(num_strata * sizeof(uint8_t *)));
    for (int i = 0; i < num_strata; i++) {
      strata[i] = static_cast<uint8_t *>(std::aligned_alloc(4096, page_aligned_size));
    }
  }

  inline ~ByteStratifiedColumn() {
    for (int i = 0; i < num_strata; i++) {
      free(strata[i]);
    }
    free(strata);
  }
};

template<class T>
void SetColumn(ByteStratifiedColumn &column, const std::vector<T> &input_data) {
  if (column.size < input_data.size()) {
    throw std::out_of_range("column size < input vec size");
  }
#pragma omp parallel for default(none) shared(column, input_data) schedule(static)
  for (size_t i = 0; i < input_data.size(); i++) {
    const T in_value = (input_data[i] << column.right_padding);
    for (int j = 0; j < column.num_strata; j++) {
      column.strata[j][i] = (in_value >> ((column.num_strata - j - 1) * 8)) & 0xFF;
    }
  }
}

template<class T, bool enable_prefetch = BYTE_FETCH_ENABLE_PREFETCH>
T SimpleFetchFromColumn(const ByteStratifiedColumn &column, size_t idx) {
  T out = 0;
  if constexpr (enable_prefetch) {
    for (int j = 0; j < column.num_strata; j++) {
      __builtin_prefetch(&column.strata[j][idx], 0, 3);
    }
  }
  for (int j = 0; j < column.num_strata; j++) {
    out |= column.strata[j][idx] << ((column.num_strata - j - 1) * 8);
  }
  out = out >> column.right_padding;
  return out;
}

template<class T, bool minimize_reads = BYTE_FETCH_MINIMIZE_READS, bool enable_prefetch = BYTE_FETCH_ENABLE_PREFETCH>
T DfeFetchFromColumn(const ByteStratifiedColumn &column, size_t idx) {
  constexpr T upper_size = enc::lib::upper_size_helper<T>();
  constexpr T upper_bytes = (upper_size / 8) + ((upper_size % 8) ? 1 : 0);
  constexpr T lower_size = (sizeof(T) * 8) - upper_size;
  T out = 0;
  if constexpr (enable_prefetch) {
    for (int j = 0; j < upper_bytes; j++) {
      __builtin_prefetch(&column.strata[j][idx], 0, 3);
    }
  }
  for (int j = 0; j < upper_bytes; j++) {
    out |= column.strata[j][idx] << ((column.num_strata - j - 1) * 8);
  }
  if (out == 0) [[unlikely]] {
    return 0;
  }
  T adtl_cols;
  if constexpr (minimize_reads) {
    adtl_cols = (out >> lower_size) - 1;
  } else {
    adtl_cols = column.num_strata - upper_size;
  }
  T total_bytes = ((upper_size + adtl_cols) / 8) + (((upper_size + adtl_cols) % 8) ? 1 : 0);
  if constexpr (enable_prefetch) {
    for (int j = upper_bytes; j < total_bytes; j++) {
      __builtin_prefetch(&column.strata[j][idx], 0, 3);
    }
  }
  for (int j = upper_bytes; j < total_bytes; j++) {
    out |= column.strata[j][idx] << ((column.num_strata - j - 1) * 8);
  }
  return enc::lib::FromDFE<T>(out);
}

template<class T, bool minimize_reads = BYTE_FETCH_MINIMIZE_READS, bool enable_prefetch = BYTE_FETCH_ENABLE_PREFETCH>
T EdfeFetchFromColumn(const ByteStratifiedColumn &column, size_t idx) {
  constexpr T upper_size = enc::lib::upper_size_helper<T>();
  constexpr T upper_bytes = ((2 + upper_size) / 8) + (((2 + upper_size) % 8) ? 1 : 0);
  constexpr T lower_size = (sizeof(T) * 8) - upper_size - 2;
  constexpr T upper_mask = (1 << upper_size) - 1;
  T out = 0;
  if constexpr (enable_prefetch) {
    for (int j = 0; j < upper_bytes; j++) {
      __builtin_prefetch(&column.strata[j][idx], 0, 3);
    }
  }
  for (int j = 0; j < upper_bytes; j++) {
    out |= column.strata[j][idx] << ((column.num_strata - j - 1) * 8);
  }
  if (out == 0) [[unlikely]] {
    return 0;
  }
  T adtl_cols;
  if constexpr (minimize_reads) {
    adtl_cols =
        ((((out >> ((sizeof(T) * 8) - 1)) ? (~out) : (out)) >> ((sizeof(T) * 8) - 2)) & 0x1)
        ? (column.num_strata - (upper_size + 2))
        : (((out >> lower_size) & (upper_mask)) - 1);
  } else {
    adtl_cols = column.num_strata - (upper_size + 2);
  }
  T total_bytes = (((2 + upper_size) + adtl_cols) / 8) + ((((2 + upper_size) + adtl_cols) % 8) ? 1 : 0);
  if constexpr (enable_prefetch) {
    for (int j = upper_bytes; j < total_bytes; j++) {
      __builtin_prefetch(&column.strata[j][idx], 0, 3);
    }
  }
  for (int j = upper_bytes; j < total_bytes; j++) {
    out |= column.strata[j][idx] << ((column.num_strata - j - 1) * 8);
  }
  return enc::lib::FromEDFE<T>(out);
}

template<
    class T,
    enc::lib::Encoding encoding,
    bool minimize_reads = BYTE_FETCH_MINIMIZE_READS,
    bool enable_prefetch = BYTE_FETCH_ENABLE_PREFETCH>
T FetchFromColumn(const ByteStratifiedColumn &column, size_t idx) {
  if constexpr (encoding == NONE) {
    return SimpleFetchFromColumn<T, enable_prefetch>(column, idx);
  }
  if constexpr (encoding == DFE) {
    return DfeFetchFromColumn<T, minimize_reads, enable_prefetch>(column, idx);
  }
  if constexpr (encoding == EDFE) {
    return EdfeFetchFromColumn<T, minimize_reads, enable_prefetch>(column, idx);
  }
}

class BitStratifiedColumn {
 public:
  const size_t size;
  const uint64_t num_strata;
  const uint64_t num_avx2_blocks;
  const uint64_t col_bytes;
  const size_t page_aligned_size;

  uint8_t **strata = nullptr;

  inline BitStratifiedColumn() = delete;
  inline BitStratifiedColumn(BitStratifiedColumn const &) = delete;
  void operator=(BitStratifiedColumn const &) = delete;

  inline BitStratifiedColumn(
      size_t column_size,
      uint64_t value_bits) : size(column_size),
                             num_strata(value_bits),
                             num_avx2_blocks((column_size / 256) + ((column_size % 256) ? 1 : 0)),
                             col_bytes((column_size / 8) + ((column_size % 8) ? 1 : 0)),
                             page_aligned_size(4096 - (col_bytes % 4096) + col_bytes) {
    strata = static_cast<uint8_t **>(std::malloc(num_strata * sizeof(uint8_t *)));
    for (int i = 0; i < num_strata; i++) {
      strata[i] = static_cast<uint8_t *>(std::aligned_alloc(4096, page_aligned_size));
    }
  }

  inline ~BitStratifiedColumn() {
    for (int i = 0; i < num_strata; i++) {
      free(strata[i]);
    }
    free(strata);
  }
};

template<class T>
void SetColumn(BitStratifiedColumn &column, const std::vector<T> &input_data) {
  if (column.size < input_data.size()) {
    throw std::out_of_range("column size < input vec size");
  }
#pragma omp parallel for default(none) shared(column, input_data) schedule(static)
  for (size_t i = 0; i < (input_data.size() / 8); i++) {
    for (int j = 0; j < column.num_strata; j++) {
      uint8_t a_val = 0;
      for (int k = 0; k < 8; k++) {
        const T extracted_bit = (input_data[(i * 8) + k] >> (column.num_strata - j - 1)) & 1;
        a_val |= extracted_bit << k;
      }
      column.strata[j][i] = a_val;
    }
  }
}

template<class T, bool enable_prefetch = BIT_FETCH_ENABLE_PREFETCH>
T SimpleFetchFromColumn(const BitStratifiedColumn &column, size_t idx) {
  T out = 0;
  if constexpr (enable_prefetch) {
    for (int j = 0; j < column.num_strata; j++) {
      __builtin_prefetch(&column.strata[j][idx / 8], 0, 3);
    }
  }
  for (int j = 0; j < column.num_strata; j++) {
    T extracted_bit = (column.strata[j][idx / 8] >> (idx % 8)) & 1;
    out |= extracted_bit << (column.num_strata - j - 1);
  }
  return out;
}

template<class T, bool minimize_reads = BIT_FETCH_MINIMIZE_READS, bool enable_prefetch = BIT_FETCH_ENABLE_PREFETCH>
T DfeFetchFromColumn(const BitStratifiedColumn &column, size_t idx) {
  constexpr T upper_size = enc::lib::upper_size_helper<T>();
  constexpr T lower_size = (sizeof(T) * 8) - upper_size;
  T out = 0;
  if constexpr (enable_prefetch) {
    for (int j = 0; j < upper_size; j++) {
      __builtin_prefetch(&column.strata[j][idx / 8], 0, 3);
    }
  }
  for (int j = 0; j < upper_size; j++) {
    T extracted_bit = (column.strata[j][idx / 8] >> (idx % 8)) & 1;
    out |= extracted_bit << (column.num_strata - j - 1);
  }
  if (out == 0) [[unlikely]] {
    return 0;
  }
  T adtl_cols;
  if constexpr (minimize_reads) {
    adtl_cols = (out >> lower_size) - 1;
  } else {
    adtl_cols = column.num_strata - upper_size;
  }
  if constexpr (enable_prefetch) {
    for (int j = upper_size; j < upper_size + adtl_cols; j++) {
      __builtin_prefetch(&column.strata[j][idx / 8], 0, 3);
    }
  }
  for (int j = upper_size; j < upper_size + adtl_cols; j++) {
    T extracted_bit = (column.strata[j][idx / 8] >> (idx % 8)) & 1;
    out |= extracted_bit << (column.num_strata - j - 1);
  }
  return enc::lib::FromDFE<T>(out);
}

template<class T, bool minimize_reads = BIT_FETCH_MINIMIZE_READS, bool enable_prefetch = BIT_FETCH_ENABLE_PREFETCH>
T EdfeFetchFromColumn(const BitStratifiedColumn &column, size_t idx) {
  constexpr T upper_size = enc::lib::upper_size_helper<T>();
  constexpr T lower_size = (sizeof(T) * 8) - upper_size - 2;
  constexpr T upper_mask = (1 << upper_size) - 1;
  T out = 0;
  if constexpr (enable_prefetch) {
    for (int j = 0; j < upper_size + 2; j++) {
      __builtin_prefetch(&column.strata[j][idx / 8], 0, 3);
    }
  }
  for (int j = 0; j < upper_size + 2; j++) {
    T extracted_bit = (column.strata[j][idx / 8] >> (idx % 8)) & 1;
    out |= extracted_bit << (column.num_strata - j - 1);
  }
  if (out == 0) [[unlikely]] {
    return 0;
  }
  T adtl_cols;
  T temp = ((out >> lower_size) & (upper_mask));
  if constexpr (minimize_reads) {
    adtl_cols =
        ((((out >> ((sizeof(T) * 8) - 1)) ? (~out) : (out)) >> ((sizeof(T) * 8) - 2)) & 0x1)
        ? (column.num_strata - (upper_size + 2))
        : (temp - 1);
  } else {
    adtl_cols = column.num_strata - (upper_size + 2);
  }
  if constexpr (enable_prefetch) {
    for (int j = upper_size + 2; j < upper_size + 2 + adtl_cols; j++) {
      __builtin_prefetch(&column.strata[j][idx / 8], 0, 3);
    }
  }
  for (int j = upper_size + 2; j < upper_size + 2 + adtl_cols; j++) {
    T extracted_bit = (column.strata[j][idx / 8] >> (idx % 8)) & 1;
    out |= extracted_bit << (column.num_strata - j - 1);
  }
  return enc::lib::FromEDFE<T>(out);
}

template<
    class T,
    enc::lib::Encoding encoding,
    bool minimize_reads = BIT_FETCH_MINIMIZE_READS,
    bool enable_prefetch = BIT_FETCH_ENABLE_PREFETCH>
T FetchFromColumn(const BitStratifiedColumn &column, size_t idx) {
  if constexpr (encoding == NONE) {
    return SimpleFetchFromColumn<T, enable_prefetch>(column, idx);
  }
  if constexpr (encoding == DFE) {
    return DfeFetchFromColumn<T, minimize_reads, enable_prefetch>(column, idx);
  }
  if constexpr (encoding == EDFE) {
    return EdfeFetchFromColumn<T, minimize_reads, enable_prefetch>(column, idx);
  }
}

enum BVModifier {
  bvSet,
  bvSetNot,
  bvAnd,
  bvAndNot,
  bvOr,
  bvOrNot
};

class Bitvector {
 public:
  const uint64_t size;
  const uint64_t blocks;
  const uint64_t page_aligned_size;
  int8_t *const bitvector;

  inline Bitvector() = delete;
  inline Bitvector(Bitvector const &) = delete;
  inline explicit Bitvector(uint64_t size) : size(size),
                                             blocks(size / 64),
                                             page_aligned_size(4096 - ((size / 8) % 4096) + (size / 8)),
                                             bitvector(static_cast<int8_t *>(std::aligned_alloc(4096, page_aligned_size))) {
    zero();
  }
  inline explicit Bitvector(const BitStratifiedColumn &column) : Bitvector(column.size){};
  inline explicit Bitvector(const ByteStratifiedColumn &column) : Bitvector(column.size){};
  inline explicit Bitvector(BitStratifiedColumn &column) : Bitvector(column.size){};
  inline explicit Bitvector(ByteStratifiedColumn &column) : Bitvector(column.size){};

  inline ~Bitvector() {
    free(bitvector);
  }

  inline void zero() {
#pragma omp parallel for schedule(static) default(none) shared(size, bitvector)
    for (int i = 0; i < (size / 64); i++) {
      __builtin_prefetch(bitvector + (i * 8) + BITVECTOR_PREFETCH_DISTANCE, 1, 3);
      *((uint64_t *) (bitvector + (i * 8))) = 0;
    }
  }

  inline uint64_t popcnt() {
    uint64_t out = 0;
#pragma omp parallel for schedule(static) default(none) shared(size, bitvector) reduction(+:out)
    for (int i = 0; i < (size / 64); i++) {
      __builtin_prefetch(bitvector + (i * 8) + BITVECTOR_PREFETCH_DISTANCE, 0, 3);
      out += _mm_popcnt_u64(*((uint64_t *) (bitvector + (i * 8))));
    }
    return out;
  }

  inline void intersect_with(const Bitvector &bv2) {
#pragma omp parallel for schedule(static) default(none) shared(size, bitvector, bv2)
    for (int i = 0; i < (size / 64); i++) {
      __builtin_prefetch(bitvector + (i * 8) + BITVECTOR_PREFETCH_DISTANCE, 1, 3);
      __builtin_prefetch(bv2.bitvector + (i * 8) + BITVECTOR_PREFETCH_DISTANCE, 0, 3);
      *((uint64_t *) (bitvector + (i * 8))) = (*((uint64_t *) (bitvector + (i * 8)))) & (*((uint64_t *) (bv2.bitvector + (i * 8))));
    }
  }
};

template<class T, enc::lib::Encoding encoding, BVModifier bv_modifier = bvSet, bool enable_prefetch = SCAN_ENABLE_PREFETCH>
inline void LTConstantForward(const ByteStratifiedColumn &column, Bitvector &result, T unpadded_constant) {
  constexpr T upper_size = lib::upper_size_helper<T>();
  constexpr T lower_size_dfe = (sizeof(T) * 8) - upper_size;
  constexpr T lower_size_edfe = (sizeof(T) * 8) - upper_size - 2;
  constexpr T upper_mask = (1 << upper_size) - 1;
  size_t end_idx;
  T constant;
  if constexpr (encoding == lib::NONE) {
    constant = unpadded_constant << column.right_padding;
    end_idx = column.num_strata - 1;
  }
  if constexpr (encoding == lib::DFE) {
    constant = unpadded_constant;
    const T end_bits = ((constant >> lower_size_dfe) & upper_mask) - 1 + upper_size;
    end_idx = ((end_bits / 8) + ((end_bits % 8) ? 1 : 0)) - 1;
  }
  if constexpr (encoding == lib::EDFE) {
    constant = unpadded_constant;
    const T end_bits = ((((constant >> ((sizeof(T) * 8) - 1)) ? (~constant) : (constant)) >> ((sizeof(T) * 8) - 2)) & 0x1)
                       ? (column.num_strata - 1)
                       : ((((constant >> lower_size_edfe) & (upper_mask)) - 1) + upper_size + 2);
    end_idx = ((end_bits / 8) + ((end_bits % 8) ? 1 : 0)) - 1;
  }
  const __m256i sign_correction = _mm256_set1_epi8((char) 0x80);
  __m256i predicate_arr[sizeof(T)];
  for (int j = 0; j <= end_idx; j++) {
    predicate_arr[j] = _mm256_xor_si256(
        sign_correction,
        _mm256_set1_epi8((char) ((constant >> ((column.num_strata - j - 1) * 8)) & 0xFF)));
  }
#pragma omp parallel for default(none) schedule(static) shared(column, result, sign_correction, predicate_arr, end_idx)
  for (int i = 0; i < column.num_avx2_blocks; i++) {
    const size_t block_idx = i * 32;
    const size_t bv_idx = i * 4;
    int eq_lock;
    if constexpr (bv_modifier == bvAnd || bv_modifier == bvAndNot) {
      eq_lock = *((int *) (result.bitvector + bv_idx));
      if (eq_lock == 0) {
        continue;
      }
    } else {
      eq_lock = -1;
    }
    int out = 0;
    for (int j = 0; j <= end_idx; j++) {
      if constexpr (enable_prefetch) {
        __builtin_prefetch(column.strata[j] + block_idx + SCAN_VERTICAL_PREFETCH_DISTANCE, 0, 3);
      }
      const __m256i data = _mm256_xor_si256(
          sign_correction,
          _mm256_lddqu_si256(reinterpret_cast<const __m256i *>(column.strata[j] + block_idx)));
      out = out | (eq_lock & _mm256_movemask_epi8(_mm256_cmpgt_epi8(predicate_arr[j], data)));
      eq_lock = eq_lock & _mm256_movemask_epi8(_mm256_cmpeq_epi8(predicate_arr[j], data));
      if (eq_lock == 0) {
        break;
      }
    }
    if constexpr (enable_prefetch) {
      __builtin_prefetch(result.bitvector + bv_idx + SCAN_VERTICAL_PREFETCH_DISTANCE, 1, 3);
    }
    if constexpr (bv_modifier == bvSet) {
      *(((int *) (result.bitvector + bv_idx))) = out;
    }
    if constexpr (bv_modifier == bvSetNot) {
      *(((int *) (result.bitvector + bv_idx))) = ~out;
    }
    if constexpr (bv_modifier == bvAnd) {
      *(((int *) (result.bitvector + bv_idx))) = *((int *) (result.bitvector + bv_idx)) & out;
    }
    if constexpr (bv_modifier == bvAndNot) {
      *(((int *) (result.bitvector + bv_idx))) = *((int *) (result.bitvector + bv_idx)) & ~out;
    }
    if constexpr (bv_modifier == bvOr) {
      *(((int *) (result.bitvector + bv_idx))) = *((int *) (result.bitvector + bv_idx)) | out;
    }
    if constexpr (bv_modifier == bvOrNot) {
      *(((int *) (result.bitvector + bv_idx))) = *((int *) (result.bitvector + bv_idx)) | ~out;
    }
  }
}

template<class T, enc::lib::Encoding encoding, BVModifier bv_modifier = bvSet, bool enable_prefetch = SCAN_ENABLE_PREFETCH>
inline void LTConstantReverse(const ByteStratifiedColumn &column, Bitvector &result, T unpadded_constant) {
  constexpr T upper_size = lib::upper_size_helper<T>();
  constexpr T lower_size_dfe = (sizeof(T) * 8) - upper_size;
  constexpr T lower_size_edfe = (sizeof(T) * 8) - upper_size - 2;
  constexpr T upper_mask = (1 << upper_size) - 1;
  size_t start_idx;
  T constant;
  if constexpr (encoding == lib::NONE) {
    constant = unpadded_constant << column.right_padding;
    start_idx = column.num_strata - 1;
  }
  if constexpr (encoding == lib::DFE) {
    constant = unpadded_constant;
    const T start_bits = ((constant >> lower_size_dfe) & upper_mask) - 1 + upper_size;
    start_idx = ((start_bits / 8) + ((start_bits % 8) ? 1 : 0)) - 1;
  }
  if constexpr (encoding == lib::EDFE) {
    constant = unpadded_constant;
    const T start_bits = ((((constant >> ((sizeof(T) * 8) - 1)) ? (~constant) : (constant)) >> ((sizeof(T) * 8) - 2)) & 0x1)
                         ? (column.num_strata - 1)
                         : ((((constant >> lower_size_edfe) & (upper_mask)) - 1) + upper_size + 2);
    start_idx = ((start_bits / 8) + ((start_bits % 8) ? 1 : 0)) - 1;
  }
  const __m256i sign_correction = _mm256_set1_epi8((char) 0x80);
  __m256i predicate_arr[sizeof(T)];
  for (int j = start_idx; j >= 0; j--) {
    predicate_arr[j] = _mm256_xor_si256(
        sign_correction,
        _mm256_set1_epi8((char) ((constant >> ((column.num_strata - j - 1) * 8)) & 0xFF)));
  }
#pragma omp parallel for default(none) schedule(static) shared(column, result, sign_correction, predicate_arr, start_idx)
  for (int i = 0; i < column.num_avx2_blocks; i++) {
    const size_t block_idx = i * 32;
    const size_t bv_idx = i * 4;
    if constexpr (bv_modifier == bvAnd || bv_modifier == bvAndNot) {
      if (*((int *) (result.bitvector + bv_idx)) == 0) {
        continue;
      }
    }
    __m256i out = _mm256_setzero_si256();
    for (int j = start_idx; j >= 0; j--) {
      if constexpr (enable_prefetch) {
        __builtin_prefetch(column.strata[j] + block_idx + SCAN_VERTICAL_PREFETCH_DISTANCE, 0, 3);
      }
      const __m256i data = _mm256_xor_si256(
          sign_correction,
          _mm256_lddqu_si256(reinterpret_cast<const __m256i *>(column.strata[j] + block_idx)));
      out = _mm256_or_si256(
          _mm256_and_si256(
              _mm256_cmpeq_epi8(predicate_arr[j], data),
              out),
          _mm256_cmpgt_epi8(predicate_arr[j], data));
    }
    if constexpr (enable_prefetch) {
      __builtin_prefetch(result.bitvector + bv_idx + SCAN_VERTICAL_PREFETCH_DISTANCE, 1, 3);
    }
    if constexpr (bv_modifier == bvSet) {
      *(((int *) (result.bitvector + bv_idx))) = _mm256_movemask_epi8(out);
    }
    if constexpr (bv_modifier == bvSetNot) {
      *(((int *) (result.bitvector + bv_idx))) = ~_mm256_movemask_epi8(out);
    }
    if constexpr (bv_modifier == bvAnd) {
      *(((int *) (result.bitvector + bv_idx))) = *((int *) (result.bitvector + bv_idx)) & _mm256_movemask_epi8(out);
    }
    if constexpr (bv_modifier == bvAndNot) {
      *(((int *) (result.bitvector + bv_idx))) = *((int *) (result.bitvector + bv_idx)) & ~_mm256_movemask_epi8(out);
    }
    if constexpr (bv_modifier == bvOr) {
      *(((int *) (result.bitvector + bv_idx))) = *((int *) (result.bitvector + bv_idx)) | _mm256_movemask_epi8(out);
    }
    if constexpr (bv_modifier == bvOrNot) {
      *(((int *) (result.bitvector + bv_idx))) = *((int *) (result.bitvector + bv_idx)) | ~_mm256_movemask_epi8(out);
    }
  }
}

template<class T, enc::lib::Encoding encoding, BVModifier bv_modifier = bvSet, bool enable_prefetch = SCAN_ENABLE_PREFETCH>
inline void LTConstantForward(const BitStratifiedColumn &column, Bitvector &result, T constant) {
  constexpr T upper_size = lib::upper_size_helper<T>();
  constexpr T lower_size_dfe = (sizeof(T) * 8) - upper_size;
  constexpr T lower_size_edfe = (sizeof(T) * 8) - upper_size - 2;
  constexpr T upper_mask = (1 << upper_size) - 1;
  size_t end_idx;
  if constexpr (encoding == lib::NONE) {
    end_idx = column.num_strata - 1;
  }
  if constexpr (encoding == lib::DFE) {
    end_idx = ((constant >> lower_size_dfe) & upper_mask) - 1 + upper_size;
  }
  if constexpr (encoding == lib::EDFE) {
    end_idx = ((((constant >> ((sizeof(T) * 8) - 1)) ? (~constant) : (constant)) >> ((sizeof(T) * 8) - 2)) & 0x1)
              ? (column.num_strata - 1)
              : ((((constant >> lower_size_edfe) & (upper_mask)) - 1) + upper_size + 2);
  }
  const __m256i ones = _mm256_set1_epi8(-1);
  __m256i predicate_arr[sizeof(T) * 8];
  for (int j = 0; j <= end_idx; j++) {
    predicate_arr[j] = (constant & (1 << (column.num_strata - j - 1))) ? ones : _mm256_setzero_si256();
  }
#pragma omp parallel for default(none) schedule(static) shared(column, result, ones, predicate_arr, end_idx)
  for (int i = 0; i < column.num_avx2_blocks; i++) {
    const size_t block_idx = i * 32;
    __m256i eq_lock;
    if constexpr (bv_modifier == bvAnd || bv_modifier == bvAndNot) {
      eq_lock = _mm256_lddqu_si256(reinterpret_cast<const __m256i *>(result.bitvector + block_idx));
      if (_mm256_testz_si256(eq_lock, eq_lock)) {
        continue;
      }
    } else {
      eq_lock = ones;
    }
    __m256i out = _mm256_setzero_si256();
    for (int j = 0; j <= end_idx; j++) {
      if constexpr (enable_prefetch) {
        __builtin_prefetch(column.strata[j] + block_idx + SCAN_VERTICAL_PREFETCH_DISTANCE, 0, 3);
      }
      const __m256i data = _mm256_lddqu_si256(reinterpret_cast<const __m256i *>(column.strata[j] + block_idx));
      out = _mm256_or_si256(
          out,
          _mm256_and_si256(
              eq_lock,
              _mm256_andnot_si256(
                  data,
                  predicate_arr[j])));
      eq_lock = _mm256_andnot_si256(
          _mm256_xor_si256(predicate_arr[j], data),
          eq_lock);
      if (_mm256_testz_si256(eq_lock, eq_lock)) {
        break;
      }
    }
    if constexpr (enable_prefetch) {
      __builtin_prefetch(result.bitvector + block_idx + SCAN_VERTICAL_PREFETCH_DISTANCE, 1, 3);
    }
    if constexpr (bv_modifier == bvSet) { _mm256_storeu_si256(reinterpret_cast<__m256i_u *>(result.bitvector + block_idx), out); }
    if constexpr (bv_modifier == bvSetNot) { _mm256_storeu_si256(reinterpret_cast<__m256i_u *>(result.bitvector + block_idx), ~out); }
    if constexpr (bv_modifier == bvAnd) {
      _mm256_storeu_si256(reinterpret_cast<__m256i_u *>(result.bitvector + block_idx),
                          _mm256_lddqu_si256(reinterpret_cast<const __m256i *>(result.bitvector + block_idx)) & out);
    }
    if constexpr (bv_modifier == bvAndNot) {
      _mm256_storeu_si256(reinterpret_cast<__m256i_u *>(result.bitvector + block_idx),
                          _mm256_lddqu_si256(reinterpret_cast<const __m256i *>(result.bitvector + block_idx)) & ~out);
    }
    if constexpr (bv_modifier == bvOr) {
      _mm256_storeu_si256(reinterpret_cast<__m256i_u *>(result.bitvector + block_idx),
                          _mm256_lddqu_si256(reinterpret_cast<const __m256i *>(result.bitvector + block_idx)) | out);
    }
    if constexpr (bv_modifier == bvOrNot) {
      _mm256_storeu_si256(reinterpret_cast<__m256i_u *>(result.bitvector + block_idx),
                          _mm256_lddqu_si256(reinterpret_cast<const __m256i *>(result.bitvector + block_idx)) | ~out);
    }
  }
}

template<class T, Encoding encoding, BVModifier bv_modifier = bvSet, bool enable_prefetch = SCAN_ENABLE_PREFETCH>
inline void LTConstantReverse(const BitStratifiedColumn &column, Bitvector &result, T constant) {
  constexpr T upper_size = lib::upper_size_helper<T>();
  constexpr T lower_size_dfe = (sizeof(T) * 8) - upper_size;
  constexpr T lower_size_edfe = (sizeof(T) * 8) - upper_size - 2;
  constexpr T upper_mask = (1 << upper_size) - 1;
  int start_idx;
  if constexpr (encoding == lib::NONE) {
    start_idx = (int) (column.num_strata - 1);
  }
  if constexpr (encoding == lib::DFE) {
    start_idx = ((constant >> lower_size_dfe) & upper_mask) - 1 + upper_size;
  }
  if constexpr (encoding == lib::EDFE) {
    start_idx = ((((constant >> ((sizeof(T) * 8) - 1)) ? (~constant) : (constant)) >> ((sizeof(T) * 8) - 2)) & 0x1)
                ? (column.num_strata - 1)
                : ((((constant >> lower_size_edfe) & (upper_mask)) - 1) + upper_size + 2);
  }
  const __m256i ones = _mm256_set1_epi8(-1);
  __m256i predicate_arr[sizeof(T) * 8];
  for (int j = 0; j < column.num_strata; j++) {
    predicate_arr[j] = (constant & (1 << (column.num_strata - j - 1))) ? ones : _mm256_setzero_si256();
  }
#pragma omp parallel for default(none) schedule(static) shared(column, result, start_idx, predicate_arr)
  for (int i = 0; i < column.num_avx2_blocks; i++) {
    const size_t block_idx = i * 32;
    if constexpr (bv_modifier == bvAnd || bv_modifier == bvAndNot) {
      const __m256i bv = _mm256_lddqu_si256(reinterpret_cast<const __m256i *>(result.bitvector + block_idx));
      if (_mm256_testz_si256(bv, bv)) {
        continue;
      }
    }
    __m256i out = _mm256_setzero_si256();
    if constexpr (enable_prefetch) {
      for (int j = start_idx; j >= 0; j--) {
        __builtin_prefetch(column.strata[j] + block_idx + SCAN_VERTICAL_PREFETCH_DISTANCE, 0, 3);
      }
    }
    for (int j = start_idx; j >= 0; j--) {
      const __m256i data = _mm256_lddqu_si256(reinterpret_cast<const __m256i *>(column.strata[j] + block_idx));
      out = _mm256_or_si256(
          _mm256_andnot_si256(
              _mm256_xor_si256(predicate_arr[j], data),
              out),
          _mm256_andnot_si256(data, predicate_arr[j]));
    }
    if constexpr (enable_prefetch) {
      __builtin_prefetch(result.bitvector + block_idx + SCAN_VERTICAL_PREFETCH_DISTANCE, 1, 3);
    }
    if constexpr (bv_modifier == bvSet) { _mm256_storeu_si256(reinterpret_cast<__m256i_u *>(result.bitvector + block_idx), out); }
    if constexpr (bv_modifier == bvSetNot) { _mm256_storeu_si256(reinterpret_cast<__m256i_u *>(result.bitvector + block_idx), ~out); }
    if constexpr (bv_modifier == bvAnd) {
      _mm256_storeu_si256(reinterpret_cast<__m256i_u *>(result.bitvector + block_idx),
                          _mm256_lddqu_si256(reinterpret_cast<const __m256i *>(result.bitvector + block_idx)) & out);
    }
    if constexpr (bv_modifier == bvAndNot) {
      _mm256_storeu_si256(reinterpret_cast<__m256i_u *>(result.bitvector + block_idx),
                          _mm256_lddqu_si256(reinterpret_cast<const __m256i *>(result.bitvector + block_idx)) & ~out);
    }
    if constexpr (bv_modifier == bvOr) {
      _mm256_storeu_si256(reinterpret_cast<__m256i_u *>(result.bitvector + block_idx),
                          _mm256_lddqu_si256(reinterpret_cast<const __m256i *>(result.bitvector + block_idx)) | out);
    }
    if constexpr (bv_modifier == bvOrNot) {
      _mm256_storeu_si256(reinterpret_cast<__m256i_u *>(result.bitvector + block_idx),
                          _mm256_lddqu_si256(reinterpret_cast<const __m256i *>(result.bitvector + block_idx)) | ~out);
    }
  }
}

template<class T, Encoding encoding, BVModifier bv_modifier = bvSet, bool forward = true, bool enable_prefetch = SCAN_ENABLE_PREFETCH>
inline void LTConstant(const BitStratifiedColumn &column, Bitvector &result, T constant) {
  if constexpr (forward) {
    LTConstantForward<T, encoding, bv_modifier, enable_prefetch>(column, result, constant);
  } else {
    LTConstantReverse<T, encoding, bv_modifier, enable_prefetch>(column, result, constant);
  }
}

template<class T, Encoding encoding, BVModifier bv_modifier = bvSet, bool forward = true, bool enable_prefetch = SCAN_ENABLE_PREFETCH>
inline void LTConstant(const ByteStratifiedColumn &column, Bitvector &result, T constant) {
  if constexpr (forward) {
    LTConstantForward<T, encoding, bv_modifier, enable_prefetch>(column, result, constant);
  } else {
    LTConstantReverse<T, encoding, bv_modifier, enable_prefetch>(column, result, constant);
  }
}

template<class T, enc::lib::Encoding encoding, BVModifier bv_modifier = bvSet, bool enable_prefetch = SCAN_ENABLE_PREFETCH>
inline void LEQConstantForward(const ByteStratifiedColumn &column, Bitvector &result, T unpadded_constant) {
  constexpr T upper_size = lib::upper_size_helper<T>();
  constexpr T lower_size_dfe = (sizeof(T) * 8) - upper_size;
  constexpr T lower_size_edfe = (sizeof(T) * 8) - upper_size - 2;
  constexpr T upper_mask = (1 << upper_size) - 1;
  size_t end_idx;
  T constant;
  if constexpr (encoding == lib::NONE) {
    constant = unpadded_constant << column.right_padding;
    end_idx = column.num_strata - 1;
  }
  if constexpr (encoding == lib::DFE) {
    constant = unpadded_constant;
    const T end_bits = ((constant >> lower_size_dfe) & upper_mask) - 1 + upper_size;
    end_idx = ((end_bits / 8) + ((end_bits % 8) ? 1 : 0)) - 1;
  }
  if constexpr (encoding == lib::EDFE) {
    constant = unpadded_constant;
    const T end_bits = ((((constant >> ((sizeof(T) * 8) - 1)) ? (~constant) : (constant)) >> ((sizeof(T) * 8) - 2)) & 0x1)
                       ? (column.num_strata - 1)
                       : ((((constant >> lower_size_edfe) & (upper_mask)) - 1) + upper_size + 2);
    end_idx = ((end_bits / 8) + ((end_bits % 8) ? 1 : 0)) - 1;
  }
  const __m256i sign_correction = _mm256_set1_epi8((char) 0x80);
  __m256i predicate_arr[sizeof(T)];
  for (int j = 0; j <= end_idx; j++) {
    predicate_arr[j] = _mm256_xor_si256(
        sign_correction,
        _mm256_set1_epi8((char) ((constant >> ((column.num_strata - j - 1) * 8)) & 0xFF)));
  }
#pragma omp parallel for default(none) schedule(static) shared(column, result, sign_correction, predicate_arr, end_idx)
  for (int i = 0; i < column.num_avx2_blocks; i++) {
    const size_t block_idx = i * 32;
    const size_t bv_idx = i * 4;
    int eq_lock;
    if constexpr (bv_modifier == bvAnd || bv_modifier == bvAndNot) {
      eq_lock = *((int *) (result.bitvector + bv_idx));
      if (eq_lock == 0) {
        continue;
      }
    } else {
      eq_lock = -1;
    }
    int out = 0;
    for (int j = 0; j <= end_idx; j++) {
      if constexpr (enable_prefetch) {
        __builtin_prefetch(column.strata[j] + block_idx + SCAN_VERTICAL_PREFETCH_DISTANCE, 0, 3);
      }
      const __m256i data = _mm256_xor_si256(
          sign_correction,
          _mm256_lddqu_si256(reinterpret_cast<const __m256i *>(column.strata[j] + block_idx)));
      out = out | (eq_lock & _mm256_movemask_epi8(_mm256_cmpgt_epi8(predicate_arr[j], data)));
      eq_lock = eq_lock & _mm256_movemask_epi8(_mm256_cmpeq_epi8(predicate_arr[j], data));
      if (eq_lock == 0) {
        break;
      }
    }
    if constexpr (enable_prefetch) {
      __builtin_prefetch(result.bitvector + bv_idx + SCAN_VERTICAL_PREFETCH_DISTANCE, 1, 3);
    }
    if constexpr (bv_modifier == bvSet) {
      *(((int *) (result.bitvector + bv_idx))) = (out | eq_lock);
    }
    if constexpr (bv_modifier == bvSetNot) {
      *(((int *) (result.bitvector + bv_idx))) = ~(out | eq_lock);
    }
    if constexpr (bv_modifier == bvAnd) {
      *(((int *) (result.bitvector + bv_idx))) = *((int *) (result.bitvector + bv_idx)) & (out | eq_lock);
    }
    if constexpr (bv_modifier == bvAndNot) {
      *(((int *) (result.bitvector + bv_idx))) = *((int *) (result.bitvector + bv_idx)) & ~(out | eq_lock);
    }
    if constexpr (bv_modifier == bvOr) {
      *(((int *) (result.bitvector + bv_idx))) = *((int *) (result.bitvector + bv_idx)) | (out | eq_lock);
    }
    if constexpr (bv_modifier == bvOrNot) {
      *(((int *) (result.bitvector + bv_idx))) = *((int *) (result.bitvector + bv_idx)) | ~(out | eq_lock);
    }
  }
}

template<class T, enc::lib::Encoding encoding, BVModifier bv_modifier = bvSet, bool enable_prefetch = SCAN_ENABLE_PREFETCH>
inline void LEQConstantReverse(const ByteStratifiedColumn &column, Bitvector &result, T unpadded_constant) {
  constexpr T upper_size = lib::upper_size_helper<T>();
  constexpr T lower_size_dfe = (sizeof(T) * 8) - upper_size;
  constexpr T lower_size_edfe = (sizeof(T) * 8) - upper_size - 2;
  constexpr T upper_mask = (1 << upper_size) - 1;
  size_t start_idx;
  T constant;
  if constexpr (encoding == lib::NONE) {
    constant = unpadded_constant << column.right_padding;
    start_idx = column.num_strata - 1;
  }
  if constexpr (encoding == lib::DFE) {
    constant = unpadded_constant;
    const T start_bits = ((constant >> lower_size_dfe) & upper_mask) - 1 + upper_size;
    start_idx = ((start_bits / 8) + ((start_bits % 8) ? 1 : 0)) - 1;
  }
  if constexpr (encoding == lib::EDFE) {
    constant = unpadded_constant;
    const T start_bits = ((((constant >> ((sizeof(T) * 8) - 1)) ? (~constant) : (constant)) >> ((sizeof(T) * 8) - 2)) & 0x1)
                         ? (column.num_strata - 1)
                         : ((((constant >> lower_size_edfe) & (upper_mask)) - 1) + upper_size + 2);
    start_idx = ((start_bits / 8) + ((start_bits % 8) ? 1 : 0)) - 1;
  }
  const __m256i sign_correction = _mm256_set1_epi8((char) 0x80);
  __m256i predicate_arr[sizeof(T)];
  for (int j = start_idx; j >= 0; j--) {
    predicate_arr[j] = _mm256_xor_si256(
        sign_correction,
        _mm256_set1_epi8((char) ((constant >> ((column.num_strata - j - 1) * 8)) & 0xFF)));
  }
#pragma omp parallel for default(none) schedule(static) shared(column, result, sign_correction, predicate_arr, start_idx)
  for (int i = 0; i < column.num_avx2_blocks; i++) {
    const size_t block_idx = i * 32;
    const size_t bv_idx = i * 4;
    if constexpr (bv_modifier == bvAnd || bv_modifier == bvAndNot) {
      if (*((int *) (result.bitvector + bv_idx)) == 0) {
        continue;
      }
    }
    int out = 0;
    int eq_lock = -1;
    for (int j = start_idx; j >= 0; j--) {
      if constexpr (enable_prefetch) {
        __builtin_prefetch(column.strata[j] + block_idx + SCAN_VERTICAL_PREFETCH_DISTANCE, 0, 3);
      }
      const __m256i data = _mm256_xor_si256(
          sign_correction,
          _mm256_lddqu_si256(reinterpret_cast<const __m256i *>(column.strata[j] + block_idx)));
      out =
          ((out) & (_mm256_movemask_epi8(_mm256_cmpeq_epi8(predicate_arr[j], data))))
              | _mm256_movemask_epi8(_mm256_cmpgt_epi8(predicate_arr[j], data));
      eq_lock =
          eq_lock & _mm256_movemask_epi8(_mm256_cmpeq_epi8(predicate_arr[j], data));
    }
    if constexpr (enable_prefetch) {
      __builtin_prefetch(result.bitvector + bv_idx + SCAN_VERTICAL_PREFETCH_DISTANCE, 1, 3);
    }
    if constexpr (bv_modifier == bvSet) {
      *(((int *) (result.bitvector + bv_idx))) = out | eq_lock;
    }
    if constexpr (bv_modifier == bvSetNot) {
      *(((int *) (result.bitvector + bv_idx))) = ~(out | eq_lock);
    }
    if constexpr (bv_modifier == bvAnd) {
      *(((int *) (result.bitvector + bv_idx))) = *((int *) (result.bitvector + bv_idx)) & (out | eq_lock);
    }
    if constexpr (bv_modifier == bvAndNot) {
      *(((int *) (result.bitvector + bv_idx))) = *((int *) (result.bitvector + bv_idx)) & ~(out | eq_lock);
    }
    if constexpr (bv_modifier == bvOr) {
      *(((int *) (result.bitvector + bv_idx))) = *((int *) (result.bitvector + bv_idx)) | (out | eq_lock);
    }
    if constexpr (bv_modifier == bvOrNot) {
      *(((int *) (result.bitvector + bv_idx))) = *((int *) (result.bitvector + bv_idx)) | ~(out | eq_lock);
    }
  }
}

template<class T, enc::lib::Encoding encoding, BVModifier bv_modifier = bvSet, bool enable_prefetch = SCAN_ENABLE_PREFETCH>
inline void LEQConstantForward(const BitStratifiedColumn &column, Bitvector &result, T constant) {
  constexpr T upper_size = lib::upper_size_helper<T>();
  constexpr T lower_size_dfe = (sizeof(T) * 8) - upper_size;
  constexpr T lower_size_edfe = (sizeof(T) * 8) - upper_size - 2;
  constexpr T upper_mask = (1 << upper_size) - 1;
  size_t end_idx;
  if constexpr (encoding == lib::NONE) {
    end_idx = column.num_strata - 1;
  }
  if constexpr (encoding == lib::DFE) {
    end_idx = ((constant >> lower_size_dfe) & upper_mask) - 1 + upper_size;
  }
  if constexpr (encoding == lib::EDFE) {
    end_idx = ((((constant >> ((sizeof(T) * 8) - 1)) ? (~constant) : (constant)) >> ((sizeof(T) * 8) - 2)) & 0x1)
              ? (column.num_strata - 1)
              : ((((constant >> lower_size_edfe) & (upper_mask)) - 1) + upper_size + 2);
  }
  const __m256i ones = _mm256_set1_epi8(-1);
  __m256i predicate_arr[sizeof(T) * 8];
  for (int j = 0; j <= end_idx; j++) {
    predicate_arr[j] = (constant & (1 << (column.num_strata - j - 1))) ? ones : _mm256_setzero_si256();
  }
#pragma omp parallel for default(none) schedule(static) shared(column, result, ones, predicate_arr, end_idx)
  for (int i = 0; i < column.num_avx2_blocks; i++) {
    const size_t block_idx = i * 32;
    __m256i eq_lock;
    if constexpr (bv_modifier == bvAnd || bv_modifier == bvAndNot) {
      eq_lock = _mm256_lddqu_si256(reinterpret_cast<const __m256i *>(result.bitvector + block_idx));
      if (_mm256_testz_si256(eq_lock, eq_lock)) {
        continue;
      }
    } else {
      eq_lock = ones;
    }
    __m256i out = _mm256_setzero_si256();
    for (int j = 0; j <= end_idx; j++) {
      if constexpr (enable_prefetch) {
        __builtin_prefetch(column.strata[j] + block_idx + SCAN_VERTICAL_PREFETCH_DISTANCE, 0, 3);
      }
      const __m256i data = _mm256_lddqu_si256(reinterpret_cast<const __m256i *>(column.strata[j] + block_idx));
      out = _mm256_or_si256(
          out,
          _mm256_and_si256(
              eq_lock,
              _mm256_andnot_si256(
                  data,
                  predicate_arr[j])));
      eq_lock = _mm256_andnot_si256(
          _mm256_xor_si256(predicate_arr[j], data),
          eq_lock);
      if (_mm256_testz_si256(eq_lock, eq_lock)) {
        break;
      }
    }
    if constexpr (enable_prefetch) {
      __builtin_prefetch(result.bitvector + block_idx + SCAN_VERTICAL_PREFETCH_DISTANCE, 1, 3);
    }
    if constexpr (bv_modifier == bvSet) {
      _mm256_storeu_si256(reinterpret_cast<__m256i_u *>(result.bitvector + block_idx), _mm256_or_si256(eq_lock, out));
    }
    if constexpr (bv_modifier == bvSetNot) {
      _mm256_storeu_si256(reinterpret_cast<__m256i_u *>(result.bitvector + block_idx), ~_mm256_or_si256(eq_lock, out));
    }
    if constexpr (bv_modifier == bvAnd) {
      _mm256_storeu_si256(reinterpret_cast<__m256i_u *>(result.bitvector + block_idx),
                          _mm256_lddqu_si256(reinterpret_cast<const __m256i *>(result.bitvector + block_idx))
                              & _mm256_or_si256(eq_lock, out));
    }
    if constexpr (bv_modifier == bvAndNot) {
      _mm256_storeu_si256(reinterpret_cast<__m256i_u *>(result.bitvector + block_idx),
                          _mm256_andnot_si256(
                              _mm256_or_si256(eq_lock, out),
                              _mm256_lddqu_si256(reinterpret_cast<const __m256i *>(result.bitvector + block_idx))));
    }
    if constexpr (bv_modifier == bvOr) {
      _mm256_storeu_si256(reinterpret_cast<__m256i_u *>(result.bitvector + block_idx),
                          _mm256_lddqu_si256(reinterpret_cast<const __m256i *>(result.bitvector + block_idx))
                              | _mm256_or_si256(eq_lock, out));
    }
    if constexpr (bv_modifier == bvOrNot) {
      _mm256_storeu_si256(reinterpret_cast<__m256i_u *>(result.bitvector + block_idx),
                          _mm256_lddqu_si256(reinterpret_cast<const __m256i *>(result.bitvector + block_idx))
                              | ~_mm256_or_si256(eq_lock, out));
    }
  }
}

template<class T, Encoding encoding, BVModifier bv_modifier = bvSet, bool enable_prefetch = SCAN_ENABLE_PREFETCH>
inline void LEQConstantReverse(const BitStratifiedColumn &column, Bitvector &result, T constant) {
  constexpr T upper_size = lib::upper_size_helper<T>();
  constexpr T lower_size_dfe = (sizeof(T) * 8) - upper_size;
  constexpr T lower_size_edfe = (sizeof(T) * 8) - upper_size - 2;
  constexpr T upper_mask = (1 << upper_size) - 1;
  int start_idx;
  if constexpr (encoding == lib::NONE) {
    start_idx = (int) (column.num_strata - 1);
  }
  if constexpr (encoding == lib::DFE) {
    start_idx = ((constant >> lower_size_dfe) & upper_mask) - 1 + upper_size;
  }
  if constexpr (encoding == lib::EDFE) {
    start_idx = ((((constant >> ((sizeof(T) * 8) - 1)) ? (~constant) : (constant)) >> ((sizeof(T) * 8) - 2)) & 0x1)
                ? (column.num_strata - 1)
                : ((((constant >> lower_size_edfe) & (upper_mask)) - 1) + upper_size + 2);
  }
  const __m256i ones = _mm256_set1_epi8(-1);
  __m256i predicate_arr[sizeof(T) * 8];
  for (int j = 0; j < column.num_strata; j++) {
    predicate_arr[j] = (constant & (1 << (column.num_strata - j - 1))) ? ones : _mm256_setzero_si256();
  }
#pragma omp parallel for default(none) schedule(static) shared(column, result, start_idx, ones, predicate_arr)
  for (int i = 0; i < column.num_avx2_blocks; i++) {
    const size_t block_idx = i * 32;
    if constexpr (bv_modifier == bvAnd || bv_modifier == bvAndNot) {
      const __m256i bv = _mm256_lddqu_si256(reinterpret_cast<const __m256i *>(result.bitvector + block_idx));
      if (_mm256_testz_si256(bv, bv)) {
        continue;
      }
    }
    __m256i out = _mm256_setzero_si256();
    __m256i eq_lock = ones;
    if constexpr (enable_prefetch) {
      for (int j = start_idx; j >= 0; j--) {
        __builtin_prefetch(column.strata[j] + block_idx + SCAN_VERTICAL_PREFETCH_DISTANCE, 0, 3);
      }
    }
    for (int j = start_idx; j >= 0; j--) {
      const __m256i data = _mm256_lddqu_si256(reinterpret_cast<const __m256i *>(column.strata[j] + block_idx));
      eq_lock = _mm256_andnot_si256(
          _mm256_xor_si256(predicate_arr[j], data),
          eq_lock);
      out = _mm256_or_si256(
          _mm256_andnot_si256(
              _mm256_xor_si256(predicate_arr[j], data),
              out),
          _mm256_andnot_si256(data, predicate_arr[j]));
    }
    if constexpr (enable_prefetch) {
      __builtin_prefetch(result.bitvector + block_idx + SCAN_VERTICAL_PREFETCH_DISTANCE, 1, 3);
    }
    if constexpr (bv_modifier == bvSet) {
      _mm256_storeu_si256(reinterpret_cast<__m256i_u *>(result.bitvector + block_idx), _mm256_or_si256(eq_lock, out));
    }
    if constexpr (bv_modifier == bvSetNot) {
      _mm256_storeu_si256(reinterpret_cast<__m256i_u *>(result.bitvector + block_idx), ~_mm256_or_si256(eq_lock, out));
    }
    if constexpr (bv_modifier == bvAnd) {
      _mm256_storeu_si256(reinterpret_cast<__m256i_u *>(result.bitvector + block_idx),
                          _mm256_lddqu_si256(reinterpret_cast<const __m256i *>(result.bitvector + block_idx))
                              & _mm256_or_si256(eq_lock, out));
    }
    if constexpr (bv_modifier == bvAndNot) {
      _mm256_storeu_si256(reinterpret_cast<__m256i_u *>(result.bitvector + block_idx),
                          _mm256_andnot_si256(
                              _mm256_or_si256(eq_lock, out),
                              _mm256_lddqu_si256(reinterpret_cast<const __m256i *>(result.bitvector + block_idx))));
    }
    if constexpr (bv_modifier == bvOr) {
      _mm256_storeu_si256(reinterpret_cast<__m256i_u *>(result.bitvector + block_idx),
                          _mm256_lddqu_si256(reinterpret_cast<const __m256i *>(result.bitvector + block_idx))
                              | _mm256_or_si256(eq_lock, out));
    }
    if constexpr (bv_modifier == bvOrNot) {
      _mm256_storeu_si256(reinterpret_cast<__m256i_u *>(result.bitvector + block_idx),
                          _mm256_lddqu_si256(reinterpret_cast<const __m256i *>(result.bitvector + block_idx))
                              | ~_mm256_or_si256(eq_lock, out));
    }
  }
}

template<class T, Encoding encoding, BVModifier bv_modifier = bvSet, bool forward = true, bool enable_prefetch = SCAN_ENABLE_PREFETCH>
inline void LEQConstant(const BitStratifiedColumn &column, Bitvector &result, T constant) {
  if constexpr (forward) {
    LEQConstantForward<T, encoding, bv_modifier, enable_prefetch>(column, result, constant);
  } else {
    LEQConstantReverse<T, encoding, bv_modifier, enable_prefetch>(column, result, constant);
  }
}

template<class T, Encoding encoding, BVModifier bv_modifier = bvSet, bool forward = true, bool enable_prefetch = SCAN_ENABLE_PREFETCH>
inline void LEQConstant(const ByteStratifiedColumn &column, Bitvector &result, T constant) {
  if constexpr (forward) {
    LEQConstantForward<T, encoding, bv_modifier, enable_prefetch>(column, result, constant);
  } else {
    LEQConstantReverse<T, encoding, bv_modifier, enable_prefetch>(column, result, constant);
  }
}

template<class T, enc::lib::Encoding encoding, BVModifier bv_modifier = bvSet, bool enable_prefetch = SCAN_ENABLE_PREFETCH>
inline void EQConstantForward(const ByteStratifiedColumn &column, Bitvector &result, T unpadded_constant) {
  constexpr T upper_size = lib::upper_size_helper<T>();
  constexpr T lower_size_dfe = (sizeof(T) * 8) - upper_size;
  constexpr T lower_size_edfe = (sizeof(T) * 8) - upper_size - 2;
  constexpr T upper_mask = (1 << upper_size) - 1;
  size_t end_idx;
  T constant;
  if constexpr (encoding == lib::NONE) {
    constant = unpadded_constant << column.right_padding;
    end_idx = column.num_strata - 1;
  }
  if constexpr (encoding == lib::DFE) {
    constant = unpadded_constant;
    const T end_bits = ((constant >> lower_size_dfe) & upper_mask) - 1 + upper_size;
    end_idx = ((end_bits / 8) + ((end_bits % 8) ? 1 : 0)) - 1;
  }
  if constexpr (encoding == lib::EDFE) {
    constant = unpadded_constant;
    const T end_bits = ((((constant >> ((sizeof(T) * 8) - 1)) ? (~constant) : (constant)) >> ((sizeof(T) * 8) - 2)) & 0x1)
                       ? (column.num_strata - 1)
                       : ((((constant >> lower_size_edfe) & (upper_mask)) - 1) + upper_size + 2);
    end_idx = ((end_bits / 8) + ((end_bits % 8) ? 1 : 0)) - 1;
  }
  const __m256i sign_correction = _mm256_set1_epi8((char) 0x80);
  __m256i predicate_arr[sizeof(T)];
  for (int j = 0; j <= end_idx; j++) {
    predicate_arr[j] = _mm256_xor_si256(
        sign_correction,
        _mm256_set1_epi8((char) ((constant >> ((column.num_strata - j - 1) * 8)) & 0xFF)));
  }
#pragma omp parallel for default(none) schedule(static) shared(column, result, sign_correction, predicate_arr, end_idx)
  for (int i = 0; i < column.num_avx2_blocks; i++) {
    const size_t block_idx = i * 32;
    const size_t bv_idx = i * 4;
    int out;
    if constexpr (bv_modifier == bvAnd || bv_modifier == bvAndNot) {
      out = *((int *) (result.bitvector + bv_idx));
      if (out == 0) {
        continue;
      }
    } else {
      out = -1;
    }
    for (int j = 0; j <= end_idx; j++) {
      if constexpr (enable_prefetch) {
        __builtin_prefetch(column.strata[j] + block_idx + SCAN_VERTICAL_PREFETCH_DISTANCE, 0, 3);
      }
      const __m256i data = _mm256_xor_si256(
          sign_correction,
          _mm256_lddqu_si256(reinterpret_cast<const __m256i *>(column.strata[j] + block_idx)));
      out = out & _mm256_movemask_epi8(_mm256_cmpeq_epi8(predicate_arr[j], data));
      if (out == 0) {
        break;
      }
    }
    if constexpr (enable_prefetch) {
      __builtin_prefetch(result.bitvector + bv_idx + SCAN_VERTICAL_PREFETCH_DISTANCE, 1, 3);
    }
    if constexpr (bv_modifier == bvSet) {
      *(((int *) (result.bitvector + bv_idx))) = out;
    }
    if constexpr (bv_modifier == bvSetNot) {
      *(((int *) (result.bitvector + bv_idx))) = ~out;
    }
    if constexpr (bv_modifier == bvAnd) {
      *(((int *) (result.bitvector + bv_idx))) = *((int *) (result.bitvector + bv_idx)) & out;
    }
    if constexpr (bv_modifier == bvAndNot) {
      *(((int *) (result.bitvector + bv_idx))) = *((int *) (result.bitvector + bv_idx)) & ~out;
    }
    if constexpr (bv_modifier == bvOr) {
      *(((int *) (result.bitvector + bv_idx))) = *((int *) (result.bitvector + bv_idx)) | out;
    }
    if constexpr (bv_modifier == bvOrNot) {
      *(((int *) (result.bitvector + bv_idx))) = *((int *) (result.bitvector + bv_idx)) | ~out;
    }
  }
}

template<class T, enc::lib::Encoding encoding, BVModifier bv_modifier = bvSet, bool enable_prefetch = SCAN_ENABLE_PREFETCH>
inline void EQConstantReverse(const ByteStratifiedColumn &column, Bitvector &result, T unpadded_constant) {
  constexpr T upper_size = lib::upper_size_helper<T>();
  constexpr T lower_size_dfe = (sizeof(T) * 8) - upper_size;
  constexpr T lower_size_edfe = (sizeof(T) * 8) - upper_size - 2;
  constexpr T upper_mask = (1 << upper_size) - 1;
  size_t start_idx;
  T constant;
  if constexpr (encoding == lib::NONE) {
    constant = unpadded_constant << column.right_padding;
    start_idx = column.num_strata - 1;
  }
  if constexpr (encoding == lib::DFE) {
    constant = unpadded_constant;
    const T start_bits = ((constant >> lower_size_dfe) & upper_mask) - 1 + upper_size;
    start_idx = ((start_bits / 8) + ((start_bits % 8) ? 1 : 0)) - 1;
  }
  if constexpr (encoding == lib::EDFE) {
    constant = unpadded_constant;
    const T start_bits = ((((constant >> ((sizeof(T) * 8) - 1)) ? (~constant) : (constant)) >> ((sizeof(T) * 8) - 2)) & 0x1)
                         ? (column.num_strata - 1)
                         : ((((constant >> lower_size_edfe) & (upper_mask)) - 1) + upper_size + 2);
    start_idx = ((start_bits / 8) + ((start_bits % 8) ? 1 : 0)) - 1;
  }
  const __m256i sign_correction = _mm256_set1_epi8((char) 0x80);
  __m256i predicate_arr[sizeof(T)];
  for (int j = start_idx; j >= 0; j--) {
    predicate_arr[j] = _mm256_xor_si256(
        sign_correction,
        _mm256_set1_epi8((char) ((constant >> ((column.num_strata - j - 1) * 8)) & 0xFF)));
  }
#pragma omp parallel for default(none) schedule(static) shared(column, result, sign_correction, predicate_arr, start_idx)
  for (int i = 0; i < column.num_avx2_blocks; i++) {
    const size_t block_idx = i * 32;
    const size_t bv_idx = i * 4;
    int out;
    if constexpr (bv_modifier == bvAnd || bv_modifier == bvAndNot) {
      out = *((int *) (result.bitvector + bv_idx));
      if (out == 0) {
        continue;
      }
    } else {
      out = -1;
    }
    for (int j = start_idx; j >= 0; j--) {
      if constexpr (enable_prefetch) {
        __builtin_prefetch(column.strata[j] + block_idx + SCAN_VERTICAL_PREFETCH_DISTANCE, 0, 3);
      }
      const __m256i data = _mm256_xor_si256(
          sign_correction,
          _mm256_lddqu_si256(reinterpret_cast<const __m256i *>(column.strata[j] + block_idx)));
      out = out & _mm256_movemask_epi8(_mm256_cmpeq_epi8(predicate_arr[j], data));
      if (out == 0) {
        break;
      }
    }
    if constexpr (enable_prefetch) {
      __builtin_prefetch(result.bitvector + bv_idx + SCAN_VERTICAL_PREFETCH_DISTANCE, 1, 3);
    }
    if constexpr (bv_modifier == bvSet) {
      *(((int *) (result.bitvector + bv_idx))) = out;
    }
    if constexpr (bv_modifier == bvSetNot) {
      *(((int *) (result.bitvector + bv_idx))) = ~out;
    }
    if constexpr (bv_modifier == bvAnd) {
      *(((int *) (result.bitvector + bv_idx))) = *((int *) (result.bitvector + bv_idx)) & out;
    }
    if constexpr (bv_modifier == bvAndNot) {
      *(((int *) (result.bitvector + bv_idx))) = *((int *) (result.bitvector + bv_idx)) & ~out;
    }
    if constexpr (bv_modifier == bvOr) {
      *(((int *) (result.bitvector + bv_idx))) = *((int *) (result.bitvector + bv_idx)) | out;
    }
    if constexpr (bv_modifier == bvOrNot) {
      *(((int *) (result.bitvector + bv_idx))) = *((int *) (result.bitvector + bv_idx)) | ~out;
    }
  }
}

template<class T, enc::lib::Encoding encoding, BVModifier bv_modifier = bvSet, bool enable_prefetch = SCAN_ENABLE_PREFETCH>
inline void EQConstantForward(const BitStratifiedColumn &column, Bitvector &result, T constant) {
  constexpr T upper_size = lib::upper_size_helper<T>();
  constexpr T lower_size_dfe = (sizeof(T) * 8) - upper_size;
  constexpr T lower_size_edfe = (sizeof(T) * 8) - upper_size - 2;
  constexpr T upper_mask = (1 << upper_size) - 1;
  size_t end_idx;
  if constexpr (encoding == lib::NONE) {
    end_idx = column.num_strata - 1;
  }
  if constexpr (encoding == lib::DFE) {
    end_idx = ((constant >> lower_size_dfe) & upper_mask) - 1 + upper_size;
  }
  if constexpr (encoding == lib::EDFE) {
    end_idx = ((((constant >> ((sizeof(T) * 8) - 1)) ? (~constant) : (constant)) >> ((sizeof(T) * 8) - 2)) & 0x1)
              ? (column.num_strata - 1)
              : ((((constant >> lower_size_edfe) & (upper_mask)) - 1) + upper_size + 2);
  }
  const __m256i ones = _mm256_set1_epi8(-1);
  __m256i predicate_arr[sizeof(T) * 8];
  for (int j = 0; j <= end_idx; j++) {
    predicate_arr[j] = (constant & (1 << (column.num_strata - j - 1))) ? ones : _mm256_setzero_si256();
  }
#pragma omp parallel for default(none) schedule(static) shared(column, result, ones, predicate_arr, end_idx)
  for (int i = 0; i < column.num_avx2_blocks; i++) {
    const size_t block_idx = i * 32;
    __m256i out;
    if constexpr (bv_modifier == bvAnd || bv_modifier == bvAndNot) {
      out = _mm256_lddqu_si256(reinterpret_cast<const __m256i *>(result.bitvector + block_idx));
      if (_mm256_testz_si256(out, out)) {
        continue;
      }
    } else {
      out = ones;
    }
    for (int j = 0; j <= end_idx; j++) {
      if constexpr (enable_prefetch) {
        __builtin_prefetch(column.strata[j] + block_idx + SCAN_VERTICAL_PREFETCH_DISTANCE, 0, 3);
      }
      const __m256i data = _mm256_lddqu_si256(reinterpret_cast<const __m256i *>(column.strata[j] + block_idx));
      out = _mm256_andnot_si256(_mm256_xor_si256(data, predicate_arr[j]), out);
      if (_mm256_testz_si256(out, out)) {
        break;
      }
    }
    if constexpr (enable_prefetch) {
      __builtin_prefetch(result.bitvector + block_idx + SCAN_VERTICAL_PREFETCH_DISTANCE, 1, 3);
    }
    if constexpr (bv_modifier == bvSet) { _mm256_storeu_si256(reinterpret_cast<__m256i_u *>(result.bitvector + block_idx), out); }
    if constexpr (bv_modifier == bvSetNot) { _mm256_storeu_si256(reinterpret_cast<__m256i_u *>(result.bitvector + block_idx), ~out); }
    if constexpr (bv_modifier == bvAnd) {
      _mm256_storeu_si256(reinterpret_cast<__m256i_u *>(result.bitvector + block_idx),
                          _mm256_lddqu_si256(reinterpret_cast<const __m256i *>(result.bitvector + block_idx)) & out);
    }
    if constexpr (bv_modifier == bvAndNot) {
      _mm256_storeu_si256(reinterpret_cast<__m256i_u *>(result.bitvector + block_idx),
                          _mm256_lddqu_si256(reinterpret_cast<const __m256i *>(result.bitvector + block_idx)) & ~out);
    }
    if constexpr (bv_modifier == bvOr) {
      _mm256_storeu_si256(reinterpret_cast<__m256i_u *>(result.bitvector + block_idx),
                          _mm256_lddqu_si256(reinterpret_cast<const __m256i *>(result.bitvector + block_idx)) | out);
    }
    if constexpr (bv_modifier == bvOrNot) {
      _mm256_storeu_si256(reinterpret_cast<__m256i_u *>(result.bitvector + block_idx),
                          _mm256_lddqu_si256(reinterpret_cast<const __m256i *>(result.bitvector + block_idx)) | ~out);
    }
  }
}

template<class T, Encoding encoding, BVModifier bv_modifier = bvSet, bool enable_prefetch = SCAN_ENABLE_PREFETCH>
inline void EQConstantReverse(const BitStratifiedColumn &column, Bitvector &result, T constant) {
  constexpr T upper_size = lib::upper_size_helper<T>();
  constexpr T lower_size_dfe = (sizeof(T) * 8) - upper_size;
  constexpr T lower_size_edfe = (sizeof(T) * 8) - upper_size - 2;
  constexpr T upper_mask = (1 << upper_size) - 1;
  int start_idx;
  if constexpr (encoding == lib::NONE) {
    start_idx = (int) (column.num_strata - 1);
  }
  if constexpr (encoding == lib::DFE) {
    start_idx = ((constant >> lower_size_dfe) & upper_mask) - 1 + upper_size;
  }
  if constexpr (encoding == lib::EDFE) {
    start_idx = ((((constant >> ((sizeof(T) * 8) - 1)) ? (~constant) : (constant)) >> ((sizeof(T) * 8) - 2)) & 0x1)
                ? (column.num_strata - 1)
                : ((((constant >> lower_size_edfe) & (upper_mask)) - 1) + upper_size + 2);
  }
  const __m256i ones = _mm256_set1_epi8(-1);
  __m256i predicate_arr[sizeof(T) * 8];
  for (int j = 0; j < column.num_strata; j++) {
    predicate_arr[j] = (constant & (1 << (column.num_strata - j - 1))) ? ones : _mm256_setzero_si256();
  }
#pragma omp parallel for default(none) schedule(static) shared(column, result, start_idx, ones, predicate_arr)
  for (int i = 0; i < column.num_avx2_blocks; i++) {
    const size_t block_idx = i * 32;
    __m256i out;
    if constexpr (bv_modifier == bvAnd || bv_modifier == bvAndNot) {
      out = _mm256_lddqu_si256(reinterpret_cast<const __m256i *>(result.bitvector + block_idx));
      if (_mm256_testz_si256(out, out)) {
        continue;
      }
    } else {
      out = ones;
    }
    if constexpr (enable_prefetch) {
      for (int j = start_idx; j >= 0; j--) {
        __builtin_prefetch(column.strata[j] + block_idx + SCAN_VERTICAL_PREFETCH_DISTANCE, 0, 3);
      }
    }
    for (int j = start_idx; j >= 0; j--) {
      const __m256i data = _mm256_lddqu_si256(reinterpret_cast<const __m256i *>(column.strata[j] + block_idx));
      out = _mm256_andnot_si256(_mm256_xor_si256(data, predicate_arr[j]), out);
      if (_mm256_testz_si256(out, out)) {
        break;
      }
    }
    if constexpr (enable_prefetch) {
      __builtin_prefetch(result.bitvector + block_idx + SCAN_VERTICAL_PREFETCH_DISTANCE, 1, 3);
    }
    if constexpr (bv_modifier == bvSet) { _mm256_storeu_si256(reinterpret_cast<__m256i_u *>(result.bitvector + block_idx), out); }
    if constexpr (bv_modifier == bvSetNot) { _mm256_storeu_si256(reinterpret_cast<__m256i_u *>(result.bitvector + block_idx), ~out); }
    if constexpr (bv_modifier == bvAnd) {
      _mm256_storeu_si256(reinterpret_cast<__m256i_u *>(result.bitvector + block_idx),
                          _mm256_lddqu_si256(reinterpret_cast<const __m256i *>(result.bitvector + block_idx)) & out);
    }
    if constexpr (bv_modifier == bvAndNot) {
      _mm256_storeu_si256(reinterpret_cast<__m256i_u *>(result.bitvector + block_idx),
                          _mm256_lddqu_si256(reinterpret_cast<const __m256i *>(result.bitvector + block_idx)) & ~out);
    }
    if constexpr (bv_modifier == bvOr) {
      _mm256_storeu_si256(reinterpret_cast<__m256i_u *>(result.bitvector + block_idx),
                          _mm256_lddqu_si256(reinterpret_cast<const __m256i *>(result.bitvector + block_idx)) | out);
    }
    if constexpr (bv_modifier == bvOrNot) {
      _mm256_storeu_si256(reinterpret_cast<__m256i_u *>(result.bitvector + block_idx),
                          _mm256_lddqu_si256(reinterpret_cast<const __m256i *>(result.bitvector + block_idx)) | ~out);
    }
  }
}

template<class T, Encoding encoding, BVModifier bv_modifier = bvSet, bool forward = true, bool enable_prefetch = SCAN_ENABLE_PREFETCH>
inline void EQConstant(const BitStratifiedColumn &column, Bitvector &result, T constant) {
  if constexpr (forward) {
    EQConstantForward<T, encoding, bv_modifier, enable_prefetch>(column, result, constant);
  } else {
    EQConstantReverse<T, encoding, bv_modifier, enable_prefetch>(column, result, constant);
  }
}

template<class T, Encoding encoding, BVModifier bv_modifier = bvSet, bool forward = true, bool enable_prefetch = SCAN_ENABLE_PREFETCH>
inline void EQConstant(const ByteStratifiedColumn &column, Bitvector &result, T constant) {
  if constexpr (forward) {
    EQConstantForward<T, encoding, bv_modifier, enable_prefetch>(column, result, constant);
  } else {
    EQConstantReverse<T, encoding, bv_modifier, enable_prefetch>(column, result, constant);
  }
}

}// namespace enc::lib

#endif//ENCODINGS_ROOT_SRC_STRATIFIED_STORAGE_H_
