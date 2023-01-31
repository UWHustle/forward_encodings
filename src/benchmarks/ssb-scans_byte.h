#ifndef ENCODINGS_ROOT_SRC_BENCH_SSB_SCANS_BYTE_H_
#define ENCODINGS_ROOT_SRC_BENCH_SSB_SCANS_BYTE_H_

#include "../encodings.h"
#include "../stratified_storage.h"
#include "../wrappers/ssb_wrapper.h"
#include "../wrappers/wrapper_utilities.h"
#include "bench_util.h"
#include <chrono>
#include <cstdlib>
#include <exception>
#include <iostream>

namespace enc::bench {

class SsbWideTableHelperByte {
 public:
  lib::ByteStratifiedColumn lo_quantity;
  lib::ByteStratifiedColumn lo_extendedprice;
  lib::ByteStratifiedColumn lo_discount;
  lib::ByteStratifiedColumn lo_revenue;
  lib::ByteStratifiedColumn lo_supplycost;
  lib::ByteStratifiedColumn p_mfgr;
  lib::ByteStratifiedColumn p_category;
  lib::ByteStratifiedColumn p_brand1;
  lib::ByteStratifiedColumn s_city;
  lib::ByteStratifiedColumn s_nation;
  lib::ByteStratifiedColumn s_region;
  lib::ByteStratifiedColumn c_city;
  lib::ByteStratifiedColumn c_nation;
  lib::ByteStratifiedColumn c_region;
  lib::ByteStratifiedColumn d_year;
  lib::ByteStratifiedColumn d_yearmonthnum;
  lib::ByteStratifiedColumn d_yearmonth;
  lib::ByteStratifiedColumn d_weeknuminyear;
  const lib::SsbWideTable &table_ref;

  inline SsbWideTableHelperByte() = delete;
  inline SsbWideTableHelperByte(SsbWideTableHelperByte const &) = delete;
  void operator=(SsbWideTableHelperByte const &) = delete;

  inline explicit SsbWideTableHelperByte(
      const lib::SsbWideTable &table, const size_t size) : lo_quantity(lib::ByteStratifiedColumn(size, 32)),
                                                           lo_extendedprice(lib::ByteStratifiedColumn(size, 32)),
                                                           lo_discount(lib::ByteStratifiedColumn(size, 32)),
                                                           lo_revenue(lib::ByteStratifiedColumn(size, 32)),
                                                           lo_supplycost(lib::ByteStratifiedColumn(size, 32)),
                                                           p_mfgr(lib::ByteStratifiedColumn(size, 32)),
                                                           p_category(lib::ByteStratifiedColumn(size, 32)),
                                                           p_brand1(lib::ByteStratifiedColumn(size, 32)),
                                                           s_city(lib::ByteStratifiedColumn(size, 32)),
                                                           s_nation(lib::ByteStratifiedColumn(size, 32)),
                                                           s_region(lib::ByteStratifiedColumn(size, 32)),
                                                           c_city(lib::ByteStratifiedColumn(size, 32)),
                                                           c_nation(lib::ByteStratifiedColumn(size, 32)),
                                                           c_region(lib::ByteStratifiedColumn(size, 32)),
                                                           d_year(lib::ByteStratifiedColumn(size, 32)),
                                                           d_yearmonthnum(lib::ByteStratifiedColumn(size, 32)),
                                                           d_yearmonth(lib::ByteStratifiedColumn(size, 32)),
                                                           d_weeknuminyear(lib::ByteStratifiedColumn(size, 32)),
                                                           table_ref(table) {}

  template<class T, lib::Encoding encoding>
  inline void load_helper(lib::SsbWideTable data_table) {
    std::vector<T> vec_lo_quantity(data_table.n());
    std::vector<T> vec_lo_extendedprice(data_table.n());
    std::vector<T> vec_lo_discount(data_table.n());
    std::vector<T> vec_lo_revenue(data_table.n());
    std::vector<T> vec_lo_supplycost(data_table.n());
    std::vector<T> vec_p_mfgr(data_table.n());
    std::vector<T> vec_p_category(data_table.n());
    std::vector<T> vec_p_brand1(data_table.n());
    std::vector<T> vec_s_city(data_table.n());
    std::vector<T> vec_s_nation(data_table.n());
    std::vector<T> vec_s_region(data_table.n());
    std::vector<T> vec_c_city(data_table.n());
    std::vector<T> vec_c_nation(data_table.n());
    std::vector<T> vec_c_region(data_table.n());
    std::vector<T> vec_d_year(data_table.n());
    std::vector<T> vec_d_yearmonthnum(data_table.n());
    std::vector<T> vec_d_yearmonth(data_table.n());
    std::vector<T> vec_d_weeknuminyear(data_table.n());
    for (size_t i = 0; i < data_table.n(); i++) {
      vec_lo_quantity[i] = predicate_helper<T, encoding>(data_table.lo_quantity[i]);
      vec_lo_extendedprice[i] = predicate_helper<T, encoding>(data_table.lo_extendedprice[i]);
      vec_lo_discount[i] = predicate_helper<T, encoding>(data_table.lo_discount[i]);
      vec_lo_revenue[i] = predicate_helper<T, encoding>(data_table.lo_revenue[i]);
      vec_lo_supplycost[i] = predicate_helper<T, encoding>(data_table.lo_supplycost[i]);
      vec_p_mfgr[i] = predicate_helper<T, encoding>(data_table.p_mfgr[i]);
      vec_p_category[i] = predicate_helper<T, encoding>(data_table.p_category[i]);
      vec_p_brand1[i] = predicate_helper<T, encoding>(data_table.p_brand1[i]);
      vec_s_city[i] = predicate_helper<T, encoding>(data_table.s_city[i]);
      vec_s_nation[i] = predicate_helper<T, encoding>(data_table.s_nation[i]);
      vec_s_region[i] = predicate_helper<T, encoding>(data_table.s_region[i]);
      vec_c_city[i] = predicate_helper<T, encoding>(data_table.c_city[i]);
      vec_c_nation[i] = predicate_helper<T, encoding>(data_table.c_nation[i]);
      vec_c_region[i] = predicate_helper<T, encoding>(data_table.c_region[i]);
      vec_d_year[i] = predicate_helper<T, encoding>(data_table.d_year[i]);
      vec_d_yearmonthnum[i] = predicate_helper<T, encoding>(data_table.d_yearmonthnum[i]);
      vec_d_yearmonth[i] = predicate_helper<T, encoding>(data_table.d_yearmonth[i]);
      vec_d_weeknuminyear[i] = predicate_helper<T, encoding>(data_table.d_weeknuminyear[i]);
    }
    lib::SetColumn(lo_quantity, vec_lo_quantity);
    lib::SetColumn(lo_extendedprice, vec_lo_extendedprice);
    lib::SetColumn(lo_discount, vec_lo_discount);
    lib::SetColumn(lo_revenue, vec_lo_revenue);
    lib::SetColumn(lo_supplycost, vec_lo_supplycost);
    lib::SetColumn(p_mfgr, vec_p_mfgr);
    lib::SetColumn(p_category, vec_p_category);
    lib::SetColumn(p_brand1, vec_p_brand1);
    lib::SetColumn(s_city, vec_s_city);
    lib::SetColumn(s_nation, vec_s_nation);
    lib::SetColumn(s_region, vec_s_region);
    lib::SetColumn(c_city, vec_c_city);
    lib::SetColumn(c_nation, vec_c_nation);
    lib::SetColumn(c_region, vec_c_region);
    lib::SetColumn(d_year, vec_d_year);
    lib::SetColumn(d_yearmonthnum, vec_d_yearmonthnum);
    lib::SetColumn(d_yearmonth, vec_d_yearmonth);
    lib::SetColumn(d_weeknuminyear, vec_d_weeknuminyear);
  }
};

template<class T, lib::Encoding encoding, bool forward_scans>
inline void ssb_scans_q11_byte(
    const lib::SsbWideTable &table,
    const SsbWideTableHelperByte &helper,
    const std::string &experiment) {
  const T predicate_eq_year = predicate_helper<T, encoding>(1993);
  const T predicate_geq_discount = predicate_helper<T, encoding>(1);
  const T predicate_leq_discount = predicate_helper<T, encoding>(3);
  const T predicate_lt_quantity = predicate_helper<T, encoding>(25);
  lib::Bitvector bv = lib::Bitvector(table.n());
  const auto q_start = std::chrono::high_resolution_clock::now();
  lib::EQConstant<T, encoding, lib::bvSet>(helper.d_year, bv, predicate_eq_year);
  lib::LEQConstant<T, encoding, lib::bvAnd>(helper.lo_discount, bv, predicate_leq_discount);
  lib::LTConstant<T, encoding, lib::bvAnd>(helper.lo_quantity, bv, predicate_lt_quantity);
  lib::LTConstant<T, encoding, lib::bvAndNot>(helper.lo_discount, bv, predicate_geq_discount);
  const auto q_end = std::chrono::high_resolution_clock::now();
  std::string popcnt = std::to_string(bv.popcnt());
  std::cout << experiment + "-" + (forward_scans ? "forward" : "reverse") + ","
      + std::to_string(std::chrono::duration_cast<std::chrono::nanoseconds>(q_end - q_start).count()) + ","
      + popcnt
            << std::endl;
}

template<class T, lib::Encoding encoding, bool forward_scans>
inline void ssb_scans_q12_byte(
    const lib::SsbWideTable &table,
    const SsbWideTableHelperByte &helper,
    const std::string &experiment) {
  const T predicate_eq_yearmonthnum = predicate_helper<T, encoding>(199401);
  const T predicate_geq_discount = predicate_helper<T, encoding>(4);
  const T predicate_leq_discount = predicate_helper<T, encoding>(6);
  const T predicate_geq_quantity = predicate_helper<T, encoding>(26);
  const T predicate_leq_quantity = predicate_helper<T, encoding>(35);
  lib::Bitvector bv = lib::Bitvector(table.n());
  const auto q_start = std::chrono::high_resolution_clock::now();
  lib::EQConstant<T, encoding, lib::bvSet>(helper.d_yearmonthnum, bv, predicate_eq_yearmonthnum);
  lib::LTConstant<T, encoding, lib::bvAndNot>(helper.lo_quantity, bv, predicate_geq_quantity);
  lib::LTConstant<T, encoding, lib::bvAndNot>(helper.lo_discount, bv, predicate_geq_discount);
  lib::LEQConstant<T, encoding, lib::bvAnd>(helper.lo_discount, bv, predicate_leq_discount);
  lib::LEQConstant<T, encoding, lib::bvAnd>(helper.lo_quantity, bv, predicate_leq_quantity);
  const auto q_end = std::chrono::high_resolution_clock::now();
  std::string popcnt = std::to_string(bv.popcnt());
  std::cout << experiment + "-" + (forward_scans ? "forward" : "reverse") + ","
      + std::to_string(std::chrono::duration_cast<std::chrono::nanoseconds>(q_end - q_start).count()) + ","
      + popcnt
            << std::endl;
}

template<class T, lib::Encoding encoding, bool forward_scans>
inline void ssb_scans_q13_byte(
    const lib::SsbWideTable &table,
    const SsbWideTableHelperByte &helper,
    const std::string &experiment) {
  const T predicate_eq_weeknuminyear = predicate_helper<T, encoding>(6);
  const T predicate_eq_year = predicate_helper<T, encoding>(1994);
  const T predicate_geq_discount = predicate_helper<T, encoding>(5);
  const T predicate_leq_discount = predicate_helper<T, encoding>(7);
  const T predicate_geq_quantity = predicate_helper<T, encoding>(36);
  const T predicate_leq_quantity = predicate_helper<T, encoding>(40);
  lib::Bitvector bv = lib::Bitvector(table.n());
  const auto q_start = std::chrono::high_resolution_clock::now();
  lib::EQConstant<T, encoding, lib::bvSet>(helper.d_weeknuminyear, bv, predicate_eq_weeknuminyear);
  lib::EQConstant<T, encoding, lib::bvAnd>(helper.d_year, bv, predicate_eq_year);
  lib::LTConstant<T, encoding, lib::bvAndNot>(helper.lo_quantity, bv, predicate_geq_quantity);
  lib::LTConstant<T, encoding, lib::bvAndNot>(helper.lo_discount, bv, predicate_geq_discount);
  lib::LEQConstant<T, encoding, lib::bvAnd>(helper.lo_discount, bv, predicate_leq_discount);
  lib::LEQConstant<T, encoding, lib::bvAnd>(helper.lo_quantity, bv, predicate_leq_quantity);
  const auto q_end = std::chrono::high_resolution_clock::now();
  std::string popcnt = std::to_string(bv.popcnt());
  std::cout << experiment + "-" + (forward_scans ? "forward" : "reverse") + ","
      + std::to_string(std::chrono::duration_cast<std::chrono::nanoseconds>(q_end - q_start).count()) + ","
      + popcnt
            << std::endl;
}

template<class T, lib::Encoding encoding, bool forward_scans>
inline void ssb_scans_q21_byte(
    const lib::SsbWideTable &table,
    const SsbWideTableHelperByte &helper,
    const std::string &experiment) {
  const T predicate_eq_category = predicate_helper<T, encoding>(table.encoder.encode_p_category("MFGR#12"));
  const T predicate_eq_region = predicate_helper<T, encoding>(table.encoder.encode_region("AMERICA"));
  lib::Bitvector bv = lib::Bitvector(table.n());
  const auto q_start = std::chrono::high_resolution_clock::now();
  lib::EQConstant<T, encoding, lib::bvSet>(helper.p_category, bv, predicate_eq_category);
  lib::EQConstant<T, encoding, lib::bvAnd>(helper.s_region, bv, predicate_eq_region);
  const auto q_end = std::chrono::high_resolution_clock::now();
  std::string popcnt = std::to_string(bv.popcnt());
  std::cout << experiment + "-" + (forward_scans ? "forward" : "reverse") + ","
      + std::to_string(std::chrono::duration_cast<std::chrono::nanoseconds>(q_end - q_start).count()) + ","
      + popcnt
            << std::endl;
}

template<class T, lib::Encoding encoding, bool forward_scans>
inline void ssb_scans_q22_byte(
    const lib::SsbWideTable &table,
    const SsbWideTableHelperByte &helper,
    const std::string &experiment) {
  const T predicate_geq_brand = predicate_helper<T, encoding>(table.encoder.encode_p_brand1("MFGR#2221"));
  const T predicate_leq_brand = predicate_helper<T, encoding>(table.encoder.encode_p_brand1("MFGR#2228"));
  const T predicate_eq_region = predicate_helper<T, encoding>(table.encoder.encode_region("ASIA"));
  lib::Bitvector bv = lib::Bitvector(table.n());
  const auto q_start = std::chrono::high_resolution_clock::now();
  lib::EQConstant<T, encoding, lib::bvSet>(helper.s_region, bv, predicate_eq_region);
  lib::LEQConstant<T, encoding, lib::bvAnd>(helper.p_brand1, bv, predicate_leq_brand);
  lib::LTConstant<T, encoding, lib::bvAndNot>(helper.p_brand1, bv, predicate_geq_brand);
  const auto q_end = std::chrono::high_resolution_clock::now();
  std::string popcnt = std::to_string(bv.popcnt());
  std::cout << experiment + "-" + (forward_scans ? "forward" : "reverse") + ","
      + std::to_string(std::chrono::duration_cast<std::chrono::nanoseconds>(q_end - q_start).count()) + ","
      + popcnt
            << std::endl;
}

template<class T, lib::Encoding encoding, bool forward_scans>
inline void ssb_scans_q23_byte(
    const lib::SsbWideTable &table,
    const SsbWideTableHelperByte &helper,
    const std::string &experiment) {
  const T predicate_eq_brand = predicate_helper<T, encoding>(table.encoder.encode_p_brand1("MFGR#2339"));
  const T predicate_eq_region = predicate_helper<T, encoding>(table.encoder.encode_region("EUROPE"));
  lib::Bitvector bv = lib::Bitvector(table.n());
  const auto q_start = std::chrono::high_resolution_clock::now();
  lib::EQConstant<T, encoding, lib::bvSet>(helper.p_brand1, bv, predicate_eq_brand);
  lib::EQConstant<T, encoding, lib::bvAnd>(helper.s_region, bv, predicate_eq_region);
  const auto q_end = std::chrono::high_resolution_clock::now();
  std::string popcnt = std::to_string(bv.popcnt());
  std::cout << experiment + "-" + (forward_scans ? "forward" : "reverse") + ","
      + std::to_string(std::chrono::duration_cast<std::chrono::nanoseconds>(q_end - q_start).count()) + ","
      + popcnt
            << std::endl;
}

template<class T, lib::Encoding encoding, bool forward_scans>
inline void ssb_scans_q31_byte(
    const lib::SsbWideTable &table,
    const SsbWideTableHelperByte &helper,
    const std::string &experiment) {
  const T predicate_eq_c_region = predicate_helper<T, encoding>(table.encoder.encode_region("ASIA"));
  const T predicate_eq_s_region = predicate_helper<T, encoding>(table.encoder.encode_region("ASIA"));
  const T predicate_geq_year = predicate_helper<T, encoding>(1992);
  const T predicate_leq_year = predicate_helper<T, encoding>(1997);
  lib::Bitvector bv = lib::Bitvector(table.n());
  const auto q_start = std::chrono::high_resolution_clock::now();
  lib::EQConstant<T, encoding, lib::bvSet>(helper.c_region, bv, predicate_eq_c_region);
  lib::EQConstant<T, encoding, lib::bvAnd>(helper.s_region, bv, predicate_eq_s_region);
  lib::LEQConstant<T, encoding, lib::bvAnd>(helper.d_year, bv, predicate_leq_year);
  lib::LTConstant<T, encoding, lib::bvAndNot>(helper.d_year, bv, predicate_geq_year);
  const auto q_end = std::chrono::high_resolution_clock::now();
  std::string popcnt = std::to_string(bv.popcnt());
  std::cout << experiment + "-" + (forward_scans ? "forward" : "reverse") + ","
      + std::to_string(std::chrono::duration_cast<std::chrono::nanoseconds>(q_end - q_start).count()) + ","
      + popcnt
            << std::endl;
}

template<class T, lib::Encoding encoding, bool forward_scans>
inline void ssb_scans_q32_byte(
    const lib::SsbWideTable &table,
    const SsbWideTableHelperByte &helper,
    const std::string &experiment) {
  const T predicate_eq_c_nation = predicate_helper<T, encoding>(table.encoder.encode_nation("UNITED STATES"));
  const T predicate_eq_s_nation = predicate_helper<T, encoding>(table.encoder.encode_nation("UNITED STATES"));
  const T predicate_geq_year = predicate_helper<T, encoding>(1992);
  const T predicate_leq_year = predicate_helper<T, encoding>(1997);
  lib::Bitvector bv = lib::Bitvector(table.n());
  const auto q_start = std::chrono::high_resolution_clock::now();
  lib::EQConstant<T, encoding, lib::bvSet>(helper.c_nation, bv, predicate_eq_c_nation);
  lib::EQConstant<T, encoding, lib::bvAnd>(helper.s_nation, bv, predicate_eq_s_nation);
  lib::LEQConstant<T, encoding, lib::bvAnd>(helper.d_year, bv, predicate_leq_year);
  lib::LTConstant<T, encoding, lib::bvAndNot>(helper.d_year, bv, predicate_geq_year);
  const auto q_end = std::chrono::high_resolution_clock::now();
  std::string popcnt = std::to_string(bv.popcnt());
  std::cout << experiment + "-" + (forward_scans ? "forward" : "reverse") + ","
      + std::to_string(std::chrono::duration_cast<std::chrono::nanoseconds>(q_end - q_start).count()) + ","
      + popcnt
            << std::endl;
}

template<class T, lib::Encoding encoding, bool forward_scans>
inline void ssb_scans_q33_byte(
    const lib::SsbWideTable &table,
    const SsbWideTableHelperByte &helper,
    const std::string &experiment) {
  const T predicate_eq_c_city1 = predicate_helper<T, encoding>(table.encoder.encode_city("UNITED KI1"));
  const T predicate_eq_c_city2 = predicate_helper<T, encoding>(table.encoder.encode_city("UNITED KI5"));
  const T predicate_eq_s_city1 = predicate_helper<T, encoding>(table.encoder.encode_city("UNITED KI1"));
  const T predicate_eq_s_city2 = predicate_helper<T, encoding>(table.encoder.encode_city("UNITED KI5"));
  const T predicate_geq_year = predicate_helper<T, encoding>(1992);
  const T predicate_leq_year = predicate_helper<T, encoding>(1997);
  lib::Bitvector bv1 = lib::Bitvector(table.n());
  lib::Bitvector bv2 = lib::Bitvector(table.n());
  const auto q_start = std::chrono::high_resolution_clock::now();
  lib::EQConstant<T, encoding, lib::bvSet>(helper.c_city, bv1, predicate_eq_c_city1);
  lib::EQConstant<T, encoding, lib::bvOr>(helper.c_city, bv1, predicate_eq_c_city2);
  lib::EQConstant<T, encoding, lib::bvSet>(helper.s_city, bv2, predicate_eq_s_city1);
  lib::EQConstant<T, encoding, lib::bvOr>(helper.s_city, bv2, predicate_eq_s_city2);
  bv1.intersect_with(bv2);
  lib::LEQConstant<T, encoding, lib::bvAnd>(helper.d_year, bv1, predicate_leq_year);
  lib::LTConstant<T, encoding, lib::bvAndNot>(helper.d_year, bv1, predicate_geq_year);
  const auto q_end = std::chrono::high_resolution_clock::now();
  std::string popcnt = std::to_string(bv1.popcnt());
  std::cout << experiment + "-" + (forward_scans ? "forward" : "reverse") + ","
      + std::to_string(std::chrono::duration_cast<std::chrono::nanoseconds>(q_end - q_start).count()) + ","
      + popcnt
            << std::endl;
}

template<class T, lib::Encoding encoding, bool forward_scans>
inline void ssb_scans_q34_byte(
    const lib::SsbWideTable &table,
    const SsbWideTableHelperByte &helper,
    const std::string &experiment) {
  const T predicate_eq_c_city1 = predicate_helper<T, encoding>(table.encoder.encode_city("UNITED KI1"));
  const T predicate_eq_c_city2 = predicate_helper<T, encoding>(table.encoder.encode_city("UNITED KI5"));
  const T predicate_eq_s_city1 = predicate_helper<T, encoding>(table.encoder.encode_city("UNITED KI1"));
  const T predicate_eq_s_city2 = predicate_helper<T, encoding>(table.encoder.encode_city("UNITED KI5"));
  const T predicate_eq_yearmonth = predicate_helper<T, encoding>(table.encoder.encode_d_yearmonth("Dec1997"));
  lib::Bitvector bv1 = lib::Bitvector(table.n());
  lib::Bitvector bv2 = lib::Bitvector(table.n());
  const auto q_start = std::chrono::high_resolution_clock::now();
  lib::EQConstant<T, encoding, lib::bvSet>(helper.c_city, bv1, predicate_eq_c_city1);
  lib::EQConstant<T, encoding, lib::bvOr>(helper.c_city, bv1, predicate_eq_c_city2);
  lib::EQConstant<T, encoding, lib::bvSet>(helper.s_city, bv2, predicate_eq_s_city1);
  lib::EQConstant<T, encoding, lib::bvOr>(helper.s_city, bv2, predicate_eq_s_city2);
  bv1.intersect_with(bv2);
  lib::EQConstant<T, encoding, lib::bvAnd>(helper.d_yearmonth, bv1, predicate_eq_yearmonth);
  const auto q_end = std::chrono::high_resolution_clock::now();
  std::string popcnt = std::to_string(bv1.popcnt());
  std::cout << experiment + "-" + (forward_scans ? "forward" : "reverse") + ","
      + std::to_string(std::chrono::duration_cast<std::chrono::nanoseconds>(q_end - q_start).count()) + ","
      + popcnt
            << std::endl;
}

template<class T, lib::Encoding encoding, bool forward_scans>
inline void ssb_scans_q41_byte(
    const lib::SsbWideTable &table,
    const SsbWideTableHelperByte &helper,
    const std::string &experiment) {
  const T predicate_eq_c_region = predicate_helper<T, encoding>(table.encoder.encode_region("AMERICA"));
  const T predicate_eq_s_region = predicate_helper<T, encoding>(table.encoder.encode_region("AMERICA"));
  const T predicate_eq_mfgr1 = predicate_helper<T, encoding>(table.encoder.encode_p_mfgr("MFGR#1"));
  const T predicate_eq_mfgr2 = predicate_helper<T, encoding>(table.encoder.encode_p_mfgr("MFGR#2"));
  lib::Bitvector bv = lib::Bitvector(table.n());
  const auto q_start = std::chrono::high_resolution_clock::now();
  lib::EQConstant<T, encoding, lib::bvSet>(helper.p_mfgr, bv, predicate_eq_mfgr1);
  lib::EQConstant<T, encoding, lib::bvOr>(helper.p_mfgr, bv, predicate_eq_mfgr2);
  lib::EQConstant<T, encoding, lib::bvAnd>(helper.c_region, bv, predicate_eq_c_region);
  lib::EQConstant<T, encoding, lib::bvAnd>(helper.s_region, bv, predicate_eq_s_region);
  const auto q_end = std::chrono::high_resolution_clock::now();
  std::string popcnt = std::to_string(bv.popcnt());
  std::cout << experiment + "-" + (forward_scans ? "forward" : "reverse") + ","
      + std::to_string(std::chrono::duration_cast<std::chrono::nanoseconds>(q_end - q_start).count()) + ","
      + popcnt
            << std::endl;
}

template<class T, lib::Encoding encoding, bool forward_scans>
inline void ssb_scans_q42_byte(
    const lib::SsbWideTable &table,
    const SsbWideTableHelperByte &helper,
    const std::string &experiment) {
  const T predicate_eq_c_region = predicate_helper<T, encoding>(table.encoder.encode_region("AMERICA"));
  const T predicate_eq_s_region = predicate_helper<T, encoding>(table.encoder.encode_region("AMERICA"));
  const T predicate_eq_mfgr1 = predicate_helper<T, encoding>(table.encoder.encode_p_mfgr("MFGR#1"));
  const T predicate_eq_mfgr2 = predicate_helper<T, encoding>(table.encoder.encode_p_mfgr("MFGR#2"));
  const T predicate_eq_year1 = predicate_helper<T, encoding>(1997);
  const T predicate_eq_year2 = predicate_helper<T, encoding>(1998);
  lib::Bitvector bv1 = lib::Bitvector(table.n());
  lib::Bitvector bv2 = lib::Bitvector(table.n());
  const auto q_start = std::chrono::high_resolution_clock::now();
  lib::EQConstant<T, encoding, lib::bvSet>(helper.p_mfgr, bv1, predicate_eq_mfgr1);
  lib::EQConstant<T, encoding, lib::bvOr>(helper.p_mfgr, bv1, predicate_eq_mfgr2);
  lib::EQConstant<T, encoding, lib::bvSet>(helper.d_year, bv2, predicate_eq_year1);
  lib::EQConstant<T, encoding, lib::bvOr>(helper.d_year, bv2, predicate_eq_year2);
  bv1.intersect_with(bv2);
  lib::EQConstant<T, encoding, lib::bvAnd>(helper.c_region, bv1, predicate_eq_c_region);
  lib::EQConstant<T, encoding, lib::bvAnd>(helper.s_region, bv1, predicate_eq_s_region);
  const auto q_end = std::chrono::high_resolution_clock::now();
  std::string popcnt = std::to_string(bv1.popcnt());
  std::cout << experiment + "-" + (forward_scans ? "forward" : "reverse") + ","
      + std::to_string(std::chrono::duration_cast<std::chrono::nanoseconds>(q_end - q_start).count()) + ","
      + popcnt
            << std::endl;
}

template<class T, lib::Encoding encoding, bool forward_scans>
inline void ssb_scans_q43_byte(
    const lib::SsbWideTable &table,
    const SsbWideTableHelperByte &helper,
    const std::string &experiment) {
  const T predicate_eq_c_region = predicate_helper<T, encoding>(table.encoder.encode_region("AMERICA"));
  const T predicate_eq_s_nation = predicate_helper<T, encoding>(table.encoder.encode_nation("UNITED STATES"));
  const T predicate_eq_category = predicate_helper<T, encoding>(table.encoder.encode_p_category("MFGR#14"));
  const T predicate_eq_year1 = predicate_helper<T, encoding>(1997);
  const T predicate_eq_year2 = predicate_helper<T, encoding>(1998);
  lib::Bitvector bv = lib::Bitvector(table.n());
  const auto q_start = std::chrono::high_resolution_clock::now();
  lib::EQConstant<T, encoding, lib::bvSet>(helper.d_year, bv, predicate_eq_year1);
  lib::EQConstant<T, encoding, lib::bvOr>(helper.d_year, bv, predicate_eq_year2);
  lib::EQConstant<T, encoding, lib::bvAnd>(helper.p_category, bv, predicate_eq_category);
  lib::EQConstant<T, encoding, lib::bvAnd>(helper.s_nation, bv, predicate_eq_s_nation);
  lib::EQConstant<T, encoding, lib::bvAnd>(helper.c_region, bv, predicate_eq_c_region);
  const auto q_end = std::chrono::high_resolution_clock::now();
  std::string popcnt = std::to_string(bv.popcnt());
  std::cout << experiment + "-" + (forward_scans ? "forward" : "reverse") + ","
      + std::to_string(std::chrono::duration_cast<std::chrono::nanoseconds>(q_end - q_start).count()) + ","
      + popcnt
            << std::endl;
}

template<class T, lib::Encoding encoding>
inline void ssb_run_queries_byte(const lib::SsbWideTable &table, const std::string &experiment, size_t repeats) {
  SsbWideTableHelperByte helper = SsbWideTableHelperByte(table, table.n());
  helper.load_helper<T, encoding>(table);
  for (size_t i = 0; i < repeats; i++) {
    ssb_scans_q11_byte<T, encoding, true>(table, helper, std::to_string(i) + ",ssb-q11-byte-" + experiment);
    ssb_scans_q11_byte<T, encoding, false>(table, helper, std::to_string(i) + ",ssb-q11-byte-" + experiment);
    ssb_scans_q12_byte<T, encoding, true>(table, helper, std::to_string(i) + ",ssb-q12-byte-" + experiment);
    ssb_scans_q12_byte<T, encoding, false>(table, helper, std::to_string(i) + ",ssb-q12-byte-" + experiment);
    ssb_scans_q13_byte<T, encoding, true>(table, helper, std::to_string(i) + ",ssb-q13-byte-" + experiment);
    ssb_scans_q13_byte<T, encoding, false>(table, helper, std::to_string(i) + ",ssb-q13-byte-" + experiment);
    ssb_scans_q21_byte<T, encoding, true>(table, helper, std::to_string(i) + ",ssb-q21-byte-" + experiment);
    ssb_scans_q21_byte<T, encoding, false>(table, helper, std::to_string(i) + ",ssb-q21-byte-" + experiment);
    ssb_scans_q22_byte<T, encoding, true>(table, helper, std::to_string(i) + ",ssb-q22-byte-" + experiment);
    ssb_scans_q22_byte<T, encoding, false>(table, helper, std::to_string(i) + ",ssb-q22-byte-" + experiment);
    ssb_scans_q23_byte<T, encoding, true>(table, helper, std::to_string(i) + ",ssb-q23-byte-" + experiment);
    ssb_scans_q23_byte<T, encoding, false>(table, helper, std::to_string(i) + ",ssb-q23-byte-" + experiment);
    ssb_scans_q31_byte<T, encoding, true>(table, helper, std::to_string(i) + ",ssb-q31-byte-" + experiment);
    ssb_scans_q31_byte<T, encoding, false>(table, helper, std::to_string(i) + ",ssb-q31-byte-" + experiment);
    ssb_scans_q32_byte<T, encoding, true>(table, helper, std::to_string(i) + ",ssb-q32-byte-" + experiment);
    ssb_scans_q32_byte<T, encoding, false>(table, helper, std::to_string(i) + ",ssb-q32-byte-" + experiment);
    ssb_scans_q33_byte<T, encoding, true>(table, helper, std::to_string(i) + ",ssb-q33-byte-" + experiment);
    ssb_scans_q33_byte<T, encoding, false>(table, helper, std::to_string(i) + ",ssb-q33-byte-" + experiment);
    ssb_scans_q34_byte<T, encoding, true>(table, helper, std::to_string(i) + ",ssb-q34-byte-" + experiment);
    ssb_scans_q34_byte<T, encoding, false>(table, helper, std::to_string(i) + ",ssb-q34-byte-" + experiment);
    ssb_scans_q41_byte<T, encoding, true>(table, helper, std::to_string(i) + ",ssb-q41-byte-" + experiment);
    ssb_scans_q41_byte<T, encoding, false>(table, helper, std::to_string(i) + ",ssb-q41-byte-" + experiment);
    ssb_scans_q42_byte<T, encoding, true>(table, helper, std::to_string(i) + ",ssb-q42-byte-" + experiment);
    ssb_scans_q42_byte<T, encoding, false>(table, helper, std::to_string(i) + ",ssb-q42-byte-" + experiment);
    ssb_scans_q43_byte<T, encoding, true>(table, helper, std::to_string(i) + ",ssb-q43-byte-" + experiment);
    ssb_scans_q43_byte<T, encoding, false>(table, helper, std::to_string(i) + ",ssb-q43-byte-" + experiment);
  }
}

inline void ssb_wrapper_byte(const lib::SsbWideTable &table, size_t repeats) {
  ssb_run_queries_byte<uint32_t, lib::NONE>(table, "natural", repeats);
  ssb_run_queries_byte<uint32_t, lib::DFE>(table, "dfe", repeats);
  ssb_run_queries_byte<int32_t, lib::EDFE>(table, "edfe", repeats);
}

}// namespace enc::bench

#endif//ENCODINGS_ROOT_SRC_BENCH_SSB_SCANS_BYTE_H_
