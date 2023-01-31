#include "tpch_wrapper.h"
#include "csv.hpp"
#include "wrapper_utilities.h"

namespace enc::lib {

TpchEncoder::TpchEncoder() {
  linestatus_dictionary_ = {
      {"O", 0},
      {"F", 1}};
  returnflag_dictionary_ = {
      {"R", 0},
      {"A", 1},
      {"N", 2}};
  shipinstruct_dictionary_ = {
      {"DELIVER IN PERSON", 0},
      {"COLLECT COD", 1},
      {"NONE", 1},
      {"TAKE BACK RETURN", 1}};
  shipmode_dictionary_ = {
      {"REG AIR", 0},
      {"AIR", 1},
      {"RAIL", 2},
      {"SHIP", 3},
      {"TRUCK", 4},
      {"MAIL", 5},
      {"FOB", 6}};
}

uint8_t TpchEncoder::encode_l_linestatus(const std::string &s) const {
  auto it = linestatus_dictionary_.find(s);
  if (it == linestatus_dictionary_.end()) {
    throw std::logic_error("could not encode linestatus from " + s);
  }
  return it->second;
}

uint8_t TpchEncoder::encode_l_returnflag(const std::string &s) const {
  auto it = returnflag_dictionary_.find(s);
  if (it == returnflag_dictionary_.end()) {
    throw std::logic_error("could not encode returnflag from " + s);
  }
  return it->second;
}

uint8_t TpchEncoder::encode_l_shipinstruct(const std::string &s) const {
  auto it = shipinstruct_dictionary_.find(s);
  if (it == shipinstruct_dictionary_.end()) {
    throw std::logic_error("could not encode shipinstruct from " + s);
  }
  return it->second;
}

uint8_t TpchEncoder::encode_l_shipmode(const std::string &s) const {
  auto it = shipmode_dictionary_.find(s);
  if (it == shipmode_dictionary_.end()) {
    throw std::logic_error("could not encode shipmode from " + s);
  }
  return it->second;
}

TpchLineItemTable load_tpch_lineitem_table(const std::string &db_path, unsigned long long scale_factor, bool use_jdc_dates) {
  TpchLineItemTable table;
  table.scale_factor = scale_factor;
  table.use_jdc_dates = use_jdc_dates;
  std::vector<std::string> column_names = {
      "l_orderkey", "l_partkey", "l_suppkey", "l_linenumber",
      "l_quantity", "l_extendedprice", "l_discount", "l_tax",
      "l_returnflag", "l_linestatus", "l_shipdate", "l_commitdate",
      "l_receiptdate", "l_shipinstruct", "l_shipmode", "l_comment"};
  csv::CSVFormat csv_format;
  csv_format.no_header().delimiter('|').column_names(column_names).quote(false);
  csv::CSVReader csv_reader(db_path, csv_format);
  for (auto &csv_row : csv_reader) {
    table.l_orderkey.push_back(csv_row["l_orderkey"].get<uint32_t>());
    table.l_partkey.push_back(csv_row["l_partkey"].get<uint32_t>());
    table.l_suppkey.push_back(csv_row["l_suppkey"].get<uint32_t>());
    table.l_linenumber.push_back(csv_row["l_linenumber"].get<uint32_t>());
    table.l_quantity.push_back(csv_row["l_quantity"].get<uint8_t>());
    table.l_extendedprice.push_back(cents_decimal_to_int(csv_row["l_extendedprice"].get<std::string>(), true));
    table.l_discount.push_back(cents_decimal_to_int(csv_row["l_discount"].get<std::string>(), true));
    table.l_tax.push_back(cents_decimal_to_int(csv_row["l_tax"].get<std::string>(), true));
    table.l_returnflag.push_back(table.encoder.encode_l_returnflag(csv_row["l_returnflag"].get<std::string>()));
    table.l_linestatus.push_back(table.encoder.encode_l_linestatus(csv_row["l_linestatus"].get<std::string>()));
    table.l_shipdate.push_back(date_string_to_int(csv_row["l_shipdate"].get<std::string>(), table.use_jdc_dates));
    table.l_commitdate.push_back(date_string_to_int(csv_row["l_commitdate"].get<std::string>(), table.use_jdc_dates));
    table.l_receiptdate.push_back(date_string_to_int(csv_row["l_receiptdate"].get<std::string>(), table.use_jdc_dates));
    table.l_shipinstruct.push_back(table.encoder.encode_l_shipinstruct(csv_row["l_shipinstruct"].get<std::string>()));
    table.l_shipmode.push_back(table.encoder.encode_l_shipmode(csv_row["l_shipmode"].get<std::string>()));
    // table.l_comment.push_back(csv_row["l_comment"].get<std::string>());
  }
  table.l_orderkey.resize(scale_factor * 6'000'000);
  table.l_partkey.resize(scale_factor * 6'000'000);
  table.l_suppkey.resize(scale_factor * 6'000'000);
  table.l_linenumber.resize(scale_factor * 6'000'000);
  table.l_quantity.resize(scale_factor * 6'000'000);
  table.l_extendedprice.resize(scale_factor * 6'000'000);
  table.l_discount.resize(scale_factor * 6'000'000);
  table.l_tax.resize(scale_factor * 6'000'000);
  table.l_returnflag.resize(scale_factor * 6'000'000);
  table.l_linestatus.resize(scale_factor * 6'000'000);
  table.l_shipdate.resize(scale_factor * 6'000'000);
  table.l_commitdate.resize(scale_factor * 6'000'000);
  table.l_receiptdate.resize(scale_factor * 6'000'000);
  table.l_shipinstruct.resize(scale_factor * 6'000'000);
  table.l_shipmode.resize(scale_factor * 6'000'000);
  // table.l_comment.resize(scale_factor * 6'000'000);
  return table;
}

}// namespace enc::lib