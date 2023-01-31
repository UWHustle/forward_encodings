#ifndef ENCODINGS_ROOT_SRC_WRAPPERS_TAXI_WRAPPER_H_
#define ENCODINGS_ROOT_SRC_WRAPPERS_TAXI_WRAPPER_H_

#include <cstdlib>
#include <string>
#include <vector>

namespace enc::lib {

struct TaxiTable {
  [[nodiscard]] size_t n() const { return trip_distance.size(); }
  std::vector<int32_t> trip_distance;
  std::vector<int32_t> total_amount;
  std::vector<int32_t> tip_amount;
  std::vector<int32_t> passenger_count;
};

TaxiTable load_taxi_table(const std::string &data_dir, size_t replications=1, bool clean_min_fare=true);

}// namespace enc::lib

#endif//ENCODINGS_ROOT_SRC_WRAPPERS_TAXI_WRAPPER_H_
