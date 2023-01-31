#ifndef ENCODINGS_ROOT_SRC_BENCH_TPCH_ETE_H_
#define ENCODINGS_ROOT_SRC_BENCH_TPCH_ETE_H_

#include "../encodings.h"
#include "../stratified_storage.h"
#include "../wrappers/tpch_wrapper.h"
#include "../wrappers/wrapper_utilities.h"
#include <chrono>
#include <cstdlib>
#include <iostream>

namespace enc::bench {

template<class T, lib::Encoding encoding, bool forward_scans>
inline void tpch_ete_q6_bit(
    const std::string &experiment,
    const lib::BitStratifiedColumn &shipdate,
    const lib::BitStratifiedColumn &discount,
    const lib::BitStratifiedColumn &quantity,
    const lib::BitStratifiedColumn &extended_price,
    const bool use_jdc_dates) {
  const T predicate_lt_shipdate = predicate_helper<T, encoding>(lib::date_string_to_int("1995-01-01", use_jdc_dates));
  const T predicate_geq_shipdate = predicate_helper<T, encoding>(lib::date_string_to_int("1994-01-01", use_jdc_dates));
  const T predicate_lt_quantity = predicate_helper<T, encoding>(24);
  const T predicate_geq_discount = predicate_helper<T, encoding>(5);
  const T predicate_leq_discount = predicate_helper<T, encoding>(7);
  lib::Bitvector bitvec(shipdate);
  const auto lt_shipdate_start = std::chrono::high_resolution_clock::now();
  lib::LTConstant<T, encoding, lib::bvSet, forward_scans>(shipdate, bitvec, predicate_lt_shipdate);
  const auto lt_shipdate_end = std::chrono::high_resolution_clock::now();
  const auto geq_discount_start = std::chrono::high_resolution_clock::now();
  lib::LTConstant<T, encoding, lib::bvAndNot, forward_scans>(discount, bitvec, predicate_geq_discount);
  const auto geq_discount_end = std::chrono::high_resolution_clock::now();
  const auto lt_quantity_start = std::chrono::high_resolution_clock::now();
  lib::LTConstant<T, encoding, lib::bvAnd, forward_scans>(quantity, bitvec, predicate_lt_quantity);
  const auto lt_quantity_end = std::chrono::high_resolution_clock::now();
  const auto geq_shipdate_start = std::chrono::high_resolution_clock::now();
  lib::LTConstant<T, encoding, lib::bvAndNot, forward_scans>(shipdate, bitvec, predicate_geq_shipdate);
  const auto geq_shipdate_end = std::chrono::high_resolution_clock::now();
  const auto leq_discount_start = std::chrono::high_resolution_clock::now();
  lib::LEQConstant<T, encoding, lib::bvAnd, forward_scans>(discount, bitvec, predicate_leq_discount);
  const auto leq_discount_end = std::chrono::high_resolution_clock::now();
  const auto agg_start = std::chrono::high_resolution_clock::now();
  double agg_out = 0;
#pragma omp parallel for schedule(static) default(none) shared(bitvec, extended_price, discount) reduction(+:agg_out)
  for (size_t i = 0; i < bitvec.blocks; i++) {
    uint64_t bv_block = *((uint64_t *) (bitvec.bitvector + (i * 8)));
    __builtin_prefetch((bitvec.bitvector + (i * 8)) + 1024);
    while (bv_block) {
      const int k = __builtin_ctzll(bv_block);
      const T a_price = lib::FetchFromColumn<T, encoding>(extended_price, (i * 64) + k);
      const T a_disc = lib::FetchFromColumn<T, encoding>(discount, (i * 64) + k);
      agg_out += a_price * (((double) a_disc) / 100);
      bv_block ^= (1LL << k);
    }
  }
  const auto agg_end = std::chrono::high_resolution_clock::now();
  const auto time_lt_shipdate =
      std::chrono::duration_cast<std::chrono::nanoseconds>(
          lt_shipdate_end - lt_shipdate_start)
          .count();
  const auto time_geq_shipdate =
      std::chrono::duration_cast<std::chrono::nanoseconds>(
          geq_shipdate_end - geq_shipdate_start)
          .count();
  const auto time_lt_quantity =
      std::chrono::duration_cast<std::chrono::nanoseconds>(
          lt_quantity_end - lt_quantity_start)
          .count();
  const auto time_leq_discount =
      std::chrono::duration_cast<std::chrono::nanoseconds>(
          leq_discount_end - leq_discount_start)
          .count();
  const auto time_geq_discount =
      std::chrono::duration_cast<std::chrono::nanoseconds>(
          geq_discount_end - geq_discount_start)
          .count();
  const auto time_agg =
      std::chrono::duration_cast<std::chrono::nanoseconds>(
          agg_end - agg_start)
          .count();
  std::string popcnt = std::to_string(bitvec.popcnt());
  std::string agg_result = std::to_string((uint64_t) (agg_out * 100));
  std::cout
      << experiment
          + std::to_string(time_lt_shipdate + time_geq_discount + time_lt_quantity + time_geq_shipdate + time_leq_discount + time_agg)
          + "," + std::to_string(time_lt_shipdate)
          + "," + std::to_string(time_geq_discount)
          + "," + std::to_string(time_lt_quantity)
          + "," + std::to_string(time_geq_shipdate)
          + "," + std::to_string(time_leq_discount)
          + "," + std::to_string(time_agg)
          + "," + popcnt + "," + agg_result
      << std::endl;
}

inline void tpch_ete_q6_bit_wrapper(const enc::lib::TpchLineItemTable &table, int repeats = 3) {
  lib::BitStratifiedColumn natural_shipdate(table.n(), 32);
  lib::BitStratifiedColumn natural_discount(table.n(), 32);
  lib::BitStratifiedColumn natural_quantity(table.n(), 32);
  lib::BitStratifiedColumn natural_extended_price(table.n(), 32);
  lib::BitStratifiedColumn dfe_shipdate(table.n(), 32);
  lib::BitStratifiedColumn dfe_discount(table.n(), 32);
  lib::BitStratifiedColumn dfe_quantity(table.n(), 32);
  lib::BitStratifiedColumn dfe_extended_price(table.n(), 32);
  lib::BitStratifiedColumn edfe_shipdate(table.n(), 32);
  lib::BitStratifiedColumn edfe_discount(table.n(), 32);
  lib::BitStratifiedColumn edfe_quantity(table.n(), 32);
  lib::BitStratifiedColumn edfe_extended_price(table.n(), 32);
  std::vector<uint32_t> vec_dfe_shipdate(table.n());
  std::vector<uint32_t> vec_dfe_discount(table.n());
  std::vector<uint32_t> vec_dfe_quantity(table.n());
  std::vector<uint32_t> vec_dfe_extended_price(table.n());
  std::vector<int32_t> vec_edfe_shipdate(table.n());
  std::vector<int32_t> vec_edfe_discount(table.n());
  std::vector<int32_t> vec_edfe_quantity(table.n());
  std::vector<int32_t> vec_edfe_extended_price(table.n());
  for (size_t i = 0; i < table.n(); i++) {
    vec_dfe_shipdate[i] = lib::ToDFE<uint32_t>(table.l_shipdate[i]);
    vec_dfe_discount[i] = lib::ToDFE<uint32_t>(table.l_discount[i]);
    vec_dfe_quantity[i] = lib::ToDFE<uint32_t>(table.l_quantity[i]);
    vec_dfe_extended_price[i] = lib::ToDFE<uint32_t>(table.l_extendedprice[i]);
    vec_edfe_shipdate[i] = lib::ToEDFE<int32_t>((int32_t) table.l_shipdate[i]);
    vec_edfe_discount[i] = lib::ToEDFE<int32_t>((int32_t) table.l_discount[i]);
    vec_edfe_quantity[i] = lib::ToEDFE<int32_t>((int32_t) table.l_quantity[i]);
    vec_edfe_extended_price[i] = lib::ToEDFE<int32_t>((int32_t) table.l_extendedprice[i]);
  }
  lib::SetColumn(natural_shipdate, table.l_shipdate);
  lib::SetColumn(natural_discount, table.l_discount);
  lib::SetColumn(natural_quantity, table.l_quantity);
  lib::SetColumn(natural_extended_price, table.l_extendedprice);
  lib::SetColumn(dfe_shipdate, vec_dfe_shipdate);
  lib::SetColumn(dfe_discount, vec_dfe_discount);
  lib::SetColumn(dfe_quantity, vec_dfe_quantity);
  lib::SetColumn(dfe_extended_price, vec_dfe_extended_price);
  lib::SetColumn(edfe_shipdate, vec_edfe_shipdate);
  lib::SetColumn(edfe_discount, vec_edfe_discount);
  lib::SetColumn(edfe_quantity, vec_edfe_quantity);
  lib::SetColumn(edfe_extended_price, vec_edfe_extended_price);
  for (int i = 0; i < repeats; i++) {
    tpch_ete_q6_bit<uint32_t, lib::NONE, true>(
        std::to_string(i) + ",tpch-q6-bit-natural-forward,",
        natural_shipdate,
        natural_discount,
        natural_quantity,
        natural_extended_price,
        table.use_jdc_dates);
  }
  for (int i = 0; i < repeats; i++) {
    tpch_ete_q6_bit<uint32_t, lib::NONE, false>(
        std::to_string(i) + ",tpch-q6-bit-natural-reverse,",
        natural_shipdate,
        natural_discount,
        natural_quantity,
        natural_extended_price,
        table.use_jdc_dates);
  }
  for (int i = 0; i < repeats; i++) {
    tpch_ete_q6_bit<uint32_t, lib::DFE, true>(
        std::to_string(i) + ",tpch-q6-bit-dfe-forward,",
        dfe_shipdate,
        dfe_discount,
        dfe_quantity,
        dfe_extended_price,
        table.use_jdc_dates);
  }
  for (int i = 0; i < repeats; i++) {
    tpch_ete_q6_bit<uint32_t, lib::DFE, false>(
        std::to_string(i) + ",tpch-q6-bit-dfe-reverse,",
        dfe_shipdate,
        dfe_discount,
        dfe_quantity,
        dfe_extended_price,
        table.use_jdc_dates);
  }
  for (int i = 0; i < repeats; i++) {
    tpch_ete_q6_bit<int32_t, lib::EDFE, true>(
        std::to_string(i) + ",tpch-q6-bit-edfe-forward,",
        edfe_shipdate,
        edfe_discount,
        edfe_quantity,
        edfe_extended_price,
        table.use_jdc_dates);
  }
  for (int i = 0; i < repeats; i++) {
    tpch_ete_q6_bit<int32_t, lib::EDFE, false>(
        std::to_string(i) + ",tpch-q6-bit-edfe-reverse,",
        edfe_shipdate,
        edfe_discount,
        edfe_quantity,
        edfe_extended_price,
        table.use_jdc_dates);
  }
}

template<class T, lib::Encoding encoding, bool forward_scans>
inline void tpch_ete_q6_byte(
    const std::string &experiment,
    const lib::ByteStratifiedColumn &shipdate,
    const lib::ByteStratifiedColumn &discount,
    const lib::ByteStratifiedColumn &quantity,
    const lib::ByteStratifiedColumn &extended_price,
    const bool use_jdc_dates) {
  const T predicate_lt_shipdate = predicate_helper<T, encoding>(lib::date_string_to_int("1995-01-01", use_jdc_dates));
  const T predicate_geq_shipdate = predicate_helper<T, encoding>(lib::date_string_to_int("1994-01-01", use_jdc_dates));
  const T predicate_lt_quantity = predicate_helper<T, encoding>(24);
  const T predicate_geq_discount = predicate_helper<T, encoding>(5);
  const T predicate_leq_discount = predicate_helper<T, encoding>(7);
  lib::Bitvector bitvec(shipdate);
  const auto lt_shipdate_start = std::chrono::high_resolution_clock::now();
  lib::LTConstant<T, encoding, lib::bvSet, forward_scans>(shipdate, bitvec, predicate_lt_shipdate);
  const auto lt_shipdate_end = std::chrono::high_resolution_clock::now();
  const auto geq_discount_start = std::chrono::high_resolution_clock::now();
  lib::LTConstant<T, encoding, lib::bvAndNot, forward_scans>(discount, bitvec, predicate_geq_discount);
  const auto geq_discount_end = std::chrono::high_resolution_clock::now();
  const auto lt_quantity_start = std::chrono::high_resolution_clock::now();
  lib::LTConstant<T, encoding, lib::bvAnd, forward_scans>(quantity, bitvec, predicate_lt_quantity);
  const auto lt_quantity_end = std::chrono::high_resolution_clock::now();
  const auto geq_shipdate_start = std::chrono::high_resolution_clock::now();
  lib::LTConstant<T, encoding, lib::bvAndNot, forward_scans>(shipdate, bitvec, predicate_geq_shipdate);
  const auto geq_shipdate_end = std::chrono::high_resolution_clock::now();
  const auto leq_discount_start = std::chrono::high_resolution_clock::now();
  lib::LEQConstant<T, encoding, lib::bvAnd, forward_scans>(discount, bitvec, predicate_leq_discount);
  const auto leq_discount_end = std::chrono::high_resolution_clock::now();
  const auto agg_start = std::chrono::high_resolution_clock::now();
  double agg_out = 0;
#pragma omp parallel for schedule(static) default(none) shared(bitvec, extended_price, discount) reduction(+:agg_out)
  for (size_t i = 0; i < bitvec.blocks; i++) {
    uint64_t bv_block = *((uint64_t *) (bitvec.bitvector + (i * 8)));
    __builtin_prefetch((bitvec.bitvector + (i * 8)) + 1024);
    while (bv_block) {
      const int k = __builtin_ctzll(bv_block);
      const T a_price = lib::FetchFromColumn<T, encoding>(extended_price, (i * 64) + k);
      const T a_disc = lib::FetchFromColumn<T, encoding>(discount, (i * 64) + k);
      agg_out += a_price * (((double) a_disc) / 100);
      bv_block ^= (1LL << k);
    }
  }
  const auto agg_end = std::chrono::high_resolution_clock::now();
  const auto time_lt_shipdate =
      std::chrono::duration_cast<std::chrono::nanoseconds>(
          lt_shipdate_end - lt_shipdate_start)
          .count();
  const auto time_geq_shipdate =
      std::chrono::duration_cast<std::chrono::nanoseconds>(
          geq_shipdate_end - geq_shipdate_start)
          .count();
  const auto time_lt_quantity =
      std::chrono::duration_cast<std::chrono::nanoseconds>(
          lt_quantity_end - lt_quantity_start)
          .count();
  const auto time_leq_discount =
      std::chrono::duration_cast<std::chrono::nanoseconds>(
          leq_discount_end - leq_discount_start)
          .count();
  const auto time_geq_discount =
      std::chrono::duration_cast<std::chrono::nanoseconds>(
          geq_discount_end - geq_discount_start)
          .count();
  const auto time_agg =
      std::chrono::duration_cast<std::chrono::nanoseconds>(
          agg_end - agg_start)
          .count();
  std::string popcnt = std::to_string(bitvec.popcnt());
  std::string agg_result = std::to_string((uint64_t) (agg_out * 100));
  std::cout
      << experiment
          + std::to_string(time_lt_shipdate + time_geq_discount + time_lt_quantity + time_geq_shipdate + time_leq_discount + time_agg)
          + "," + std::to_string(time_lt_shipdate)
          + "," + std::to_string(time_geq_discount)
          + "," + std::to_string(time_lt_quantity)
          + "," + std::to_string(time_geq_shipdate)
          + "," + std::to_string(time_leq_discount)
          + "," + std::to_string(time_agg)
          + "," + popcnt + "," + agg_result
      << std::endl;
}

inline void tpch_ete_q6_byte_wrapper(const enc::lib::TpchLineItemTable &table, int repeats = 3) {
  lib::ByteStratifiedColumn natural_shipdate(table.n(), 32);
  lib::ByteStratifiedColumn natural_discount(table.n(), 32);
  lib::ByteStratifiedColumn natural_quantity(table.n(), 32);
  lib::ByteStratifiedColumn natural_extended_price(table.n(), 32);
  lib::ByteStratifiedColumn dfe_shipdate(table.n(), 32);
  lib::ByteStratifiedColumn dfe_discount(table.n(), 32);
  lib::ByteStratifiedColumn dfe_quantity(table.n(), 32);
  lib::ByteStratifiedColumn dfe_extended_price(table.n(), 32);
  lib::ByteStratifiedColumn edfe_shipdate(table.n(), 32);
  lib::ByteStratifiedColumn edfe_discount(table.n(), 32);
  lib::ByteStratifiedColumn edfe_quantity(table.n(), 32);
  lib::ByteStratifiedColumn edfe_extended_price(table.n(), 32);
  std::vector<uint32_t> vec_dfe_shipdate(table.n());
  std::vector<uint32_t> vec_dfe_discount(table.n());
  std::vector<uint32_t> vec_dfe_quantity(table.n());
  std::vector<uint32_t> vec_dfe_extended_price(table.n());
  std::vector<int32_t> vec_edfe_shipdate(table.n());
  std::vector<int32_t> vec_edfe_discount(table.n());
  std::vector<int32_t> vec_edfe_quantity(table.n());
  std::vector<int32_t> vec_edfe_extended_price(table.n());
  for (size_t i = 0; i < table.n(); i++) {
    vec_dfe_shipdate[i] = lib::ToDFE<uint32_t>(table.l_shipdate[i]);
    vec_dfe_discount[i] = lib::ToDFE<uint32_t>(table.l_discount[i]);
    vec_dfe_quantity[i] = lib::ToDFE<uint32_t>(table.l_quantity[i]);
    vec_dfe_extended_price[i] = lib::ToDFE<uint32_t>(table.l_extendedprice[i]);
    vec_edfe_shipdate[i] = lib::ToEDFE<int32_t>((int32_t) table.l_shipdate[i]);
    vec_edfe_discount[i] = lib::ToEDFE<int32_t>((int32_t) table.l_discount[i]);
    vec_edfe_quantity[i] = lib::ToEDFE<int32_t>((int32_t) table.l_quantity[i]);
    vec_edfe_extended_price[i] = lib::ToEDFE<int32_t>((int32_t) table.l_extendedprice[i]);
  }
  lib::SetColumn(natural_shipdate, table.l_shipdate);
  lib::SetColumn(natural_discount, table.l_discount);
  lib::SetColumn(natural_quantity, table.l_quantity);
  lib::SetColumn(natural_extended_price, table.l_extendedprice);
  lib::SetColumn(dfe_shipdate, vec_dfe_shipdate);
  lib::SetColumn(dfe_discount, vec_dfe_discount);
  lib::SetColumn(dfe_quantity, vec_dfe_quantity);
  lib::SetColumn(dfe_extended_price, vec_dfe_extended_price);
  lib::SetColumn(edfe_shipdate, vec_edfe_shipdate);
  lib::SetColumn(edfe_discount, vec_edfe_discount);
  lib::SetColumn(edfe_quantity, vec_edfe_quantity);
  lib::SetColumn(edfe_extended_price, vec_edfe_extended_price);
  for (int i = 0; i < repeats; i++) {
    tpch_ete_q6_byte<uint32_t, lib::NONE, true>(
        std::to_string(i) + ",tpch-q6-byte-natural-forward,",
        natural_shipdate,
        natural_discount,
        natural_quantity,
        natural_extended_price,
        table.use_jdc_dates);
  }
  for (int i = 0; i < repeats; i++) {
    tpch_ete_q6_byte<uint32_t, lib::NONE, false>(
        std::to_string(i) + ",tpch-q6-byte-natural-reverse,",
        natural_shipdate,
        natural_discount,
        natural_quantity,
        natural_extended_price,
        table.use_jdc_dates);
  }
  for (int i = 0; i < repeats; i++) {
    tpch_ete_q6_byte<uint32_t, lib::DFE, true>(
        std::to_string(i) + ",tpch-q6-byte-dfe-forward,",
        dfe_shipdate,
        dfe_discount,
        dfe_quantity,
        dfe_extended_price,
        table.use_jdc_dates);
  }
  for (int i = 0; i < repeats; i++) {
    tpch_ete_q6_byte<uint32_t, lib::DFE, false>(
        std::to_string(i) + ",tpch-q6-byte-dfe-reverse,",
        dfe_shipdate,
        dfe_discount,
        dfe_quantity,
        dfe_extended_price,
        table.use_jdc_dates);
  }
  for (int i = 0; i < repeats; i++) {
    tpch_ete_q6_byte<int32_t, lib::EDFE, true>(
        std::to_string(i) + ",tpch-q6-byte-edfe-forward,",
        edfe_shipdate,
        edfe_discount,
        edfe_quantity,
        edfe_extended_price,
        table.use_jdc_dates);
  }
  for (int i = 0; i < repeats; i++) {
    tpch_ete_q6_byte<int32_t, lib::EDFE, false>(
        std::to_string(i) + ",tpch-q6-byte-edfe-reverse,",
        edfe_shipdate,
        edfe_discount,
        edfe_quantity,
        edfe_extended_price,
        table.use_jdc_dates);
  }
}

}// namespace enc::bench

#endif//ENCODINGS_ROOT_SRC_BENCH_TPCH_ETE_H_
