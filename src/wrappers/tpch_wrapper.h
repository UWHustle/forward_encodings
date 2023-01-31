#ifndef ENCODINGS_ROOT_SRC_WRAPPERS_TPCH_WRAPPER_H_
#define ENCODINGS_ROOT_SRC_WRAPPERS_TPCH_WRAPPER_H_

#include <cstdlib>
#include <string>
#include <unordered_map>
#include <vector>

namespace enc::lib {

class TpchEncoder {
 public:
  TpchEncoder();
  [[nodiscard]] uint8_t encode_l_returnflag(const std::string &s) const;
  [[nodiscard]] uint8_t encode_l_linestatus(const std::string &s) const;
  [[nodiscard]] uint8_t encode_l_shipinstruct(const std::string &s) const;
  [[nodiscard]] uint8_t encode_l_shipmode(const std::string &s) const;

 private:
  std::unordered_map<std::string, uint8_t> returnflag_dictionary_;
  std::unordered_map<std::string, uint8_t> linestatus_dictionary_;
  std::unordered_map<std::string, uint8_t> shipinstruct_dictionary_;
  std::unordered_map<std::string, uint8_t> shipmode_dictionary_;
};

struct TpchLineItemTable {

  [[nodiscard]] size_t n() const { return l_orderkey.size(); }
  unsigned long long scale_factor;
  bool use_jdc_dates;

  std::vector<uint32_t> l_orderkey;
  std::vector<uint32_t> l_partkey;
  std::vector<uint32_t> l_suppkey;
  std::vector<uint32_t> l_linenumber;
  std::vector<uint8_t> l_quantity;
  std::vector<uint32_t> l_extendedprice;
  std::vector<uint8_t> l_discount;
  std::vector<uint8_t> l_tax;
  std::vector<uint8_t> l_returnflag;
  std::vector<uint8_t> l_linestatus;
  std::vector<uint32_t> l_shipdate;
  std::vector<uint32_t> l_commitdate;
  std::vector<uint32_t> l_receiptdate;
  std::vector<uint8_t> l_shipinstruct;
  std::vector<uint8_t> l_shipmode;
  std::vector<std::string> l_comment;

  TpchEncoder encoder;
};

TpchLineItemTable load_tpch_lineitem_table(const std::string &db_path, unsigned long long scale_factor, bool use_jdc_dates);

}// namespace enc::lib

#endif//ENCODINGS_ROOT_SRC_WRAPPERS_TPCH_WRAPPER_H_
