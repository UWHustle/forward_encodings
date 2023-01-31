#ifndef ENCODINGS_ROOT_SRC_EXPERIMENTS_COMPRESSION_H_
#define ENCODINGS_ROOT_SRC_EXPERIMENTS_COMPRESSION_H_

#include "../encodings.h"
#include "../stratified_storage.h"
#include "../wrappers/taxi_wrapper.h"
#include <cstdlib>
#include <zlib.h>

namespace enc::exp {

#define COMPRESSION_LEVEL Z_BEST_COMPRESSION

struct CompressionReport {
  size_t uncompressed = 0;
  size_t bit_natural = 0;
  size_t bit_dfe = 0;
  size_t bit_edfe = 0;
  size_t byte_natural = 0;
  size_t byte_dfe = 0;
  size_t byte_edfe = 0;
  size_t word_natural = 0;
  size_t word_dfe = 0;
  size_t word_edfe = 0;
};

template<class signed_T, class unsigned_T, class input_T>
class CompressionData {
 public:
  bool _use_dfe;
  bool _use_edfe;
  lib::BitStratifiedColumn _col_bit_natural;
  lib::BitStratifiedColumn _col_bit_dfe;
  lib::BitStratifiedColumn _col_bit_edfe;
  lib::ByteStratifiedColumn _col_byte_natural;
  lib::ByteStratifiedColumn _col_byte_dfe;
  lib::ByteStratifiedColumn _col_byte_edfe;
  const size_t _dataset_size;
  std::vector<signed_T> _dataset_natural;
  std::vector<unsigned_T> _dataset_dfe;
  std::vector<signed_T> _dataset_edfe;

  inline CompressionData() = delete;
  inline CompressionData(CompressionData const &) = delete;

  inline CompressionData(
      std::vector<input_T> dataset,
      int32_t bits_natural,
      int32_t bits_dfe,
      int32_t bits_edfe,
      bool use_dfe = true,
      bool use_edfe = true) : _col_bit_natural(lib::BitStratifiedColumn(dataset.size(), bits_natural)),
                              _col_bit_dfe(lib::BitStratifiedColumn(dataset.size(), bits_dfe)),
                              _col_bit_edfe(lib::BitStratifiedColumn(dataset.size(), bits_edfe)),
                              _col_byte_natural(lib::ByteStratifiedColumn(dataset.size(), bits_natural)),
                              _col_byte_dfe(lib::ByteStratifiedColumn(dataset.size(), bits_dfe)),
                              _col_byte_edfe(lib::ByteStratifiedColumn(dataset.size(), bits_edfe)),
                              _use_dfe(use_dfe),
                              _use_edfe(use_edfe),
                              _dataset_size(dataset.size()),
                              _dataset_natural(dataset.size()),
                              _dataset_dfe(dataset.size()),
                              _dataset_edfe(dataset.size()) {
#pragma omp parallel for schedule(static) default(none) shared(dataset, _dataset_dfe, _dataset_edfe, _use_dfe, _use_edfe)
    for (size_t i = 0; i < _dataset_size; i++) {
      _dataset_natural[i] = dataset[i];
      if (_use_dfe) {
        try {
          _dataset_dfe[i] = lib::ToDFE<unsigned_T>((unsigned_T) dataset[i]);
        } catch (...) {
          _use_dfe = false;
        }
      }
      if (_use_edfe) {
        try {
          _dataset_edfe[i] = lib::ToEDFE<signed_T>(dataset[i]);
        } catch (...) {
          _use_edfe = false;
        }
      }
    }
    lib::SetColumn(_col_bit_natural, dataset);
    lib::SetColumn(_col_byte_natural, dataset);
    if (use_dfe) {
      lib::SetColumn(_col_bit_dfe, _dataset_dfe);
      lib::SetColumn(_col_byte_dfe, _dataset_dfe);
    }
    if (use_edfe) {
      lib::SetColumn(_col_bit_edfe, _dataset_edfe);
      lib::SetColumn(_col_byte_edfe, _dataset_edfe);
    }
  }
};

template<class T>
inline size_t test_compression(const std::vector<T> &col) {
  const size_t flattened_size = col.size() * sizeof(T);
  const size_t buffer_size = compressBound(flattened_size);
  auto flattened = (T*) malloc(col.size() * sizeof(T));
#pragma omp parallel for default(none) shared(col, flattened) schedule(static)
  for(size_t i = 0; i < col.size(); i++) {
    flattened[i] = col[i];
  }
  auto *buffer = (Bytef *) malloc(buffer_size);
  uLongf dst_size = buffer_size;
  compress2(buffer, &dst_size, (Bytef *)flattened, flattened_size, COMPRESSION_LEVEL);
  free(buffer);
  free(flattened);
  return dst_size;
}

inline size_t test_compression(const lib::BitStratifiedColumn &col) {
  const size_t flattened_size = col.col_bytes * col.num_strata;
  const size_t buffer_size = compressBound(flattened_size);
  auto flattened = (uint8_t*) malloc(col.col_bytes * col.num_strata);
  for (size_t i = 0; i < col.num_strata; i++) {
#pragma omp parallel for default(none) shared(i, col, flattened) schedule(static)
    for(size_t j = 0; j < col.col_bytes; j++) {
      flattened[(i * col.col_bytes) + j] = col.strata[i][j];
    }
  }
  auto *buffer = (Bytef *) malloc(buffer_size);
  uLongf dst_size = buffer_size;
  compress2(buffer, &dst_size, (Bytef *)flattened, flattened_size, COMPRESSION_LEVEL);
  free(buffer);
  free(flattened);
  return dst_size;
}

inline size_t test_compression(const lib::ByteStratifiedColumn &col) {
  const size_t flattened_size = col.size * col.num_strata;
  const size_t buffer_size = compressBound(flattened_size);
  auto flattened = (uint8_t*) malloc(col.size * col.num_strata);
  for (size_t i = 0; i < col.num_strata; i++) {
#pragma omp parallel for default(none) shared(i, col, flattened) schedule(static)
    for(size_t j = 0; j < col.size; j++) {
      flattened[(i * col.size) + j] = col.strata[i][j];
    }
  }
  auto *buffer = (Bytef *) malloc(buffer_size);
  uLongf dst_size = buffer_size;
  compress2(buffer, &dst_size, (Bytef *)flattened, flattened_size, COMPRESSION_LEVEL);
  free(buffer);
  free(flattened);
  return dst_size;
}

template<class T>
inline void compress_column(const std::vector<T> &col, const std::string &col_name) {
  const auto data = CompressionData<int32_t, uint32_t, T>(col, 32, 32, 32);
  CompressionReport report;
  report.uncompressed = col.size() * sizeof(T);
  report.bit_natural = test_compression(data._col_bit_natural);
  report.bit_dfe = test_compression(data._col_bit_dfe);
  report.bit_edfe = test_compression(data._col_bit_edfe);
  report.byte_natural = test_compression(data._col_byte_natural);
  report.byte_dfe = test_compression(data._col_byte_dfe);
  report.byte_edfe = test_compression(data._col_byte_edfe);
  report.word_natural = test_compression(data._dataset_natural);
  report.word_dfe = test_compression(data._dataset_dfe);
  report.word_edfe = test_compression(data._dataset_edfe);
  std::cout << col_name
          + "," + std::to_string(report.uncompressed) + ","
          + std::to_string(report.bit_natural) + ","
          + std::to_string(report.bit_dfe) + ","
          + std::to_string(report.bit_edfe) + ","
          + std::to_string(report.byte_natural) + ","
          + std::to_string(report.byte_dfe) + ","
          + std::to_string(report.byte_edfe) + ","
          + std::to_string(report.word_natural) + ","
          + std::to_string(report.word_dfe) + ","
          + std::to_string(report.word_edfe)
            << std::endl;
}

inline void compression_experiment(const lib::TaxiTable &taxi) {
  std::cout << "experiment,uncompressed,"
               "bit_natural,bit_dfe,bit_edfe,"
               "byte_natural,byte_dfe,byte_edfe,"
               "word_natural,word_dfe,word_edfe"
            << std::endl;
  compress_column(taxi.passenger_count, "passenger_count");
  compress_column(taxi.tip_amount, "tip_amount");
  compress_column(taxi.trip_distance, "trip_distance");
}

}// namespace enc::exp

#endif//ENCODINGS_ROOT_SRC_EXPERIMENTS_COMPRESSION_H_
