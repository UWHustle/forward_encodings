#include "taxi_wrapper.h"
#include "wrapper_utilities.h"
#include <csv.hpp>

namespace enc::lib {

TaxiTable load_taxi_table(const std::string &db_path, const size_t replications, bool clean_min_fare) {
  TaxiTable table;
  csv::CSVFormat csv_format;
  csv::CSVReader csv_reader(db_path, csv_format);
  for (auto &csv_row : csv_reader) {
    std::string raw_total_amount = csv_row["total_amount"].get<std::string>();
    std::string raw_passenger_count = csv_row["passenger_count"].get<std::string>();
    std::string raw_trip_distance = csv_row["trip_distance"].get<std::string>();
    std::string raw_tip_amount = csv_row["tip_amount"].get<std::string>();
    // minimum fare changed from $2.50 to $3.00 on december 19th, 2022
    // https://www.nyc.gov/site/tlc/passengers/taxi-fare.page
    if (raw_total_amount.length()
        && raw_passenger_count.length()
        && raw_trip_distance.length()
        && raw_tip_amount.length()) {
      const auto total_amount = cents_decimal_to_int(raw_total_amount, true);
      if ((!clean_min_fare) || (total_amount >= 250)) {
        table.total_amount.push_back(total_amount);
        table.passenger_count.push_back(cents_decimal_to_int(raw_passenger_count, false));
        table.trip_distance.push_back(cents_decimal_to_int(raw_trip_distance, true));
        table.tip_amount.push_back(cents_decimal_to_int(raw_tip_amount, true));
      }
    }
  }
  const size_t pre_replicate_size = table.total_amount.size();
  for (size_t i = 0; i < replications; i++) {
    for (size_t j = 0; j < pre_replicate_size; j++) {
      table.total_amount.push_back(table.total_amount[j]);
      table.trip_distance.push_back(table.trip_distance[j]);
      table.passenger_count.push_back(table.passenger_count[j]);
      table.tip_amount.push_back(table.tip_amount[j]);
    }
  }
  return table;
}

}// namespace enc::lib