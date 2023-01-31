#include "ssb_wrapper.h"
#include "csv.hpp"
#include <stdexcept>

namespace enc::lib {

SsbEncoder::SsbEncoder() {
  for (uint8_t i = 0; i < 25; ++i) {
    category_dictionary_.emplace("MFGR#" + std::to_string(i / 5 + 1) + std::to_string(i % 5 + 1),
                                 i);
  }

  for (uint16_t i = 0; i < 1000; ++i) {
    brand1_dictionary_.emplace("MFGR#" + std::to_string(i / 200 + 1) + std::to_string((i % 200) / 40 + 1) + std::to_string(i % 40 + 1),
                               i);
  }

  std::vector<std::string> city_prefixes = {
      "ALGERIA  ", "ARGENTINA", "BRAZIL   ", "CANADA   ", "EGYPT    ", "ETHIOPIA ", "FRANCE   ",
      "GERMANY  ", "INDIA    ", "INDONESIA", "IRAN     ", "IRAQ     ", "JAPAN    ", "JORDAN   ",
      "KENYA    ", "MOROCCO  ", "MOZAMBIQU", "PERU     ", "CHINA    ", "ROMANIA  ", "SAUDI ARA",
      "VIETNAM  ", "RUSSIA   ", "UNITED KI", "UNITED ST"};
  for (uint8_t i = 0; i < 250; ++i) {
    city_dictionary_.emplace(city_prefixes[i / 10] + std::to_string(i % 10), i);
  }

  month_dictionary_ = {{"Jan", 0}, {"Feb", 1}, {"Mar", 2}, {"Apr", 3}, {"May", 4}, {"Jun", 5}, {"Jul", 6}, {"Aug", 7}, {"Sep", 8}, {"Oct", 9}, {"Nov", 10}, {"Dec", 11}};

  mfgr_dictionary_ = {{"MFGR#1", 0}, {"MFGR#2", 1}, {"MFGR#3", 2}, {"MFGR#4", 3}, {"MFGR#5", 4}};

  nation_dictionary_ = {
      {"ALGERIA", 0},
      {"ARGENTINA", 1},
      {"BRAZIL", 2},
      {"CANADA", 3},
      {"EGYPT", 4},
      {"ETHIOPIA", 5},
      {"FRANCE", 6},
      {"GERMANY", 7},
      {"INDIA", 8},
      {"INDONESIA", 9},
      {"IRAN", 10},
      {"IRAQ", 11},
      {"JAPAN", 12},
      {"JORDAN", 13},
      {"KENYA", 14},
      {"MOROCCO", 15},
      {"MOZAMBIQUE", 16},
      {"PERU", 17},
      {"CHINA", 18},
      {"ROMANIA", 19},
      {"SAUDI ARABIA", 20},
      {"VIETNAM", 21},
      {"RUSSIA", 22},
      {"UNITED KINGDOM", 23},
      {"UNITED STATES", 24}};

  region_dictionary_ = {
      {"AFRICA", 0},
      {"AMERICA", 1},
      {"ASIA", 2},
      {"EUROPE", 3},
      {"MIDDLE EAST", 4}};
}

uint8_t SsbEncoder::encode_d_yearmonth(const std::string &s) const {
  auto it = month_dictionary_.find(s.substr(0, 3));
  if (it == month_dictionary_.end()) {
    throw std::logic_error("could not encode yearmonth from " + s);
  }
  return (it->second << 3) | (std::stoul(s.substr(3, 4)) - 1992);
}

uint8_t SsbEncoder::encode_p_mfgr(const std::string &s) const {
  auto it = mfgr_dictionary_.find(s);
  if (it == mfgr_dictionary_.end()) {
    throw std::logic_error("could not encode p_mfgr from " + s);
  }
  return it->second;
}

uint8_t SsbEncoder::encode_p_category(const std::string &s) const {
  auto it = category_dictionary_.find(s);
  if (it == category_dictionary_.end()) {
    throw std::logic_error("could not encode p_category from " + s);
  }
  return it->second;
}

uint16_t SsbEncoder::encode_p_brand1(const std::string &s) const {
  auto it = brand1_dictionary_.find(s);
  if (it == brand1_dictionary_.end()) {
    throw std::logic_error("could not encode p_brand1 from " + s);
  }
  return it->second;
}

uint8_t SsbEncoder::encode_city(const std::string &s) const {
  auto it = city_dictionary_.find(s);
  if (it == city_dictionary_.end()) {
    throw std::logic_error("could not encode city from " + s);
  }
  return it->second;
}

uint8_t SsbEncoder::encode_nation(const std::string &s) const {
  auto it = nation_dictionary_.find(s);
  if (it == nation_dictionary_.end()) {
    throw std::logic_error("could not encode nation from " + s);
  }
  return it->second;
}

uint8_t SsbEncoder::encode_region(const std::string &s) const {
  auto it = region_dictionary_.find(s);
  if (it == region_dictionary_.end()) {
    throw std::logic_error("could not encode nation from " + s);
  }
  return it->second;
}

SsbWideTable load_ssb_widetable(const std::string &db_path, unsigned long long scale_factor) {
  SsbWideTable table;
  table.scale_factor = scale_factor;
  std::vector<std::string> column_names = {
      "lo_quantity", "lo_extendedprice", "lo_discount",
      "lo_revenue", "lo_supplycost", "p_mfgr", "p_category", "p_brand1",
      "s_city", "s_nation", "s_region", "c_city", "c_nation",
      "c_region", "d_year", "d_yearmonthnum", "d_yearmonth", "d_weeknuminyear"};
  csv::CSVFormat csv_format;
  csv_format.no_header().delimiter('|').column_names(column_names).quote(false);
  csv::CSVReader csv_reader(db_path, csv_format);
  for (auto &csv_row : csv_reader) {
    table.lo_quantity.push_back(csv_row["lo_quantity"].get<uint8_t>());
    table.lo_extendedprice.push_back(csv_row["lo_extendedprice"].get<uint32_t>());
    table.lo_discount.push_back(csv_row["lo_discount"].get<uint8_t>());
    table.lo_revenue.push_back(csv_row["lo_revenue"].get<uint32_t>());
    table.lo_supplycost.push_back(csv_row["lo_supplycost"].get<uint32_t>());
    table.p_mfgr.push_back(table.encoder.encode_p_mfgr(csv_row["p_mfgr"].get<std::string>()));
    table.p_category.push_back(table.encoder.encode_p_category(csv_row["p_category"].get<std::string>()));
    table.p_brand1.push_back(table.encoder.encode_p_brand1(csv_row["p_brand1"].get<std::string>()));
    table.s_city.push_back(table.encoder.encode_city(csv_row["s_city"].get<std::string>()));
    table.s_nation.push_back(table.encoder.encode_nation(csv_row["s_nation"].get<std::string>()));
    table.s_region.push_back(table.encoder.encode_region(csv_row["s_region"].get<std::string>()));
    table.c_city.push_back(table.encoder.encode_city(csv_row["c_city"].get<std::string>()));
    table.c_nation.push_back(table.encoder.encode_nation(csv_row["c_nation"].get<std::string>()));
    table.c_region.push_back(table.encoder.encode_region(csv_row["c_region"].get<std::string>()));
    table.d_year.push_back(csv_row["d_year"].get<uint16_t>());
    table.d_yearmonthnum.push_back(csv_row["d_yearmonthnum"].get<uint32_t>());
    table.d_yearmonth.push_back(table.encoder.encode_d_yearmonth(csv_row["d_yearmonth"].get<std::string>()));
    table.d_weeknuminyear.push_back(csv_row["d_weeknuminyear"].get<uint8_t>());
  }
  table.lo_quantity.resize(scale_factor * 6'000'000);
  table.lo_extendedprice.resize(scale_factor * 6'000'000);
  table.lo_discount.resize(scale_factor * 6'000'000);
  table.lo_revenue.resize(scale_factor * 6'000'000);
  table.lo_supplycost.resize(scale_factor * 6'000'000);
  table.p_mfgr.resize(scale_factor * 6'000'000);
  table.p_category.resize(scale_factor * 6'000'000);
  table.p_brand1.resize(scale_factor * 6'000'000);
  table.s_city.resize(scale_factor * 6'000'000);
  table.s_nation.resize(scale_factor * 6'000'000);
  table.s_region.resize(scale_factor * 6'000'000);
  table.c_city.resize(scale_factor * 6'000'000);
  table.c_nation.resize(scale_factor * 6'000'000);
  table.c_region.resize(scale_factor * 6'000'000);
  table.d_year.resize(scale_factor * 6'000'000);
  table.d_yearmonthnum.resize(scale_factor * 6'000'000);
  table.d_yearmonth.resize(scale_factor * 6'000'000);
  table.d_weeknuminyear.resize(scale_factor * 6'000'000);
  return table;
}

}// namespace enc::lib