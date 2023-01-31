#ifndef ENCODINGS_ROOT_SRC_SSB_WRAPPER_H_
#define ENCODINGS_ROOT_SRC_SSB_WRAPPER_H_

#include <cstdint>
#include <cstdlib>
#include <string>
#include <unordered_map>
#include <vector>

namespace enc::lib {

class SsbEncoder {
 public:
  SsbEncoder();
  [[nodiscard]] uint8_t encode_d_yearmonth(const std::string &s) const;
  [[nodiscard]] uint8_t encode_p_mfgr(const std::string &s) const;
  [[nodiscard]] uint8_t encode_p_category(const std::string &s) const;
  [[nodiscard]] uint16_t encode_p_brand1(const std::string &s) const;
  [[nodiscard]] uint8_t encode_city(const std::string &s) const;
  [[nodiscard]] uint8_t encode_nation(const std::string &s) const;
  [[nodiscard]] uint8_t encode_region(const std::string &s) const;

 private:
  std::unordered_map<std::string, uint8_t> category_dictionary_;
  std::unordered_map<std::string, uint16_t> brand1_dictionary_;
  std::unordered_map<std::string, uint8_t> city_dictionary_;
  std::unordered_map<std::string, uint8_t> month_dictionary_;
  std::unordered_map<std::string, uint8_t> mfgr_dictionary_;
  std::unordered_map<std::string, uint8_t> nation_dictionary_;
  std::unordered_map<std::string, uint8_t> region_dictionary_;
};

struct SsbWideTable {

  [[nodiscard]] size_t n() const { return lo_quantity.size(); }
  unsigned long long scale_factor;

  std::vector<uint8_t> lo_quantity;
  std::vector<uint32_t> lo_extendedprice;
  std::vector<uint8_t> lo_discount;
  std::vector<uint32_t> lo_revenue;
  std::vector<uint32_t> lo_supplycost;
  std::vector<uint8_t> p_mfgr;
  std::vector<uint8_t> p_category;
  std::vector<uint16_t> p_brand1;
  std::vector<uint8_t> s_city;
  std::vector<uint8_t> s_nation;
  std::vector<uint8_t> s_region;
  std::vector<uint8_t> c_city;
  std::vector<uint8_t> c_nation;
  std::vector<uint8_t> c_region;
  std::vector<uint16_t> d_year;
  std::vector<uint32_t> d_yearmonthnum;
  std::vector<uint8_t> d_yearmonth;
  std::vector<uint8_t> d_weeknuminyear;

  SsbEncoder encoder;
};

SsbWideTable load_ssb_widetable(const std::string &db_path, unsigned long long scale_factor);

}// namespace enc::lib

#endif//ENCODINGS_ROOT_SRC_SSB_WRAPPER_H_
