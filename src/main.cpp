#include "benchmarks/microbenchmarks.h"
#include "benchmarks/ssb-scans_bit.h"
#include "benchmarks/ssb-scans_byte.h"
#include "benchmarks/tpch-ete.h"
#include "experiments/compression.h"
#include "wrappers/ssb_wrapper.h"
#include "wrappers/taxi_wrapper.h"

#include <cxxopts.hpp>
#include <iostream>

int main(int argc, char **argv) {
  cxxopts::Options options("Encodings Benchmark", "DFE/EDFE Benchmark (VLDB 2023)");
  options.add_options()
      ("a,repeats", "Experiment Repeats", cxxopts::value<int>()->default_value("1"))
      ("f,db_scale_factor", "TPC-H/SSB Scale Factor", cxxopts::value<int>()->default_value("1"))
      ("r,taxi_scale_factor", "Taxi Scale Factor", cxxopts::value<int>()->default_value("1"))
      ("t,tpch_path", "TPC-H File Path", cxxopts::value<std::string>()->default_value("UNSET"))
      ("x,taxi_path", "Jan21 Yellow Trip Data File Path", cxxopts::value<std::string>()->default_value("UNSET"))
      ("w,ssb_wide_path", "SSB Wide Table Path", cxxopts::value<std::string>()->default_value("UNSET"))
      ("e,experiment", "Run compression experiment instead of benchmarks", cxxopts::value<bool>()->default_value("false")->implicit_value("true"));
  auto result = options.parse(argc, argv);
  if (result.count("help")) {
    std::cout << options.help() << std::endl;
    return 0;
  }
  int repeats = result["repeats"].as<int>();
  int db_scale_factor = result["db_scale_factor"].as<int>();
  int taxi_scale_factor = result["taxi_scale_factor"].as<int>();
  bool experiment = result["experiment"].as<bool>();
  std::string tpch_path = result["tpch_path"].as<std::string>();
  std::string ssb_wide_path = result["ssb_wide_path"].as<std::string>();
  std::string taxi_path = result["taxi_path"].as<std::string>();
  if (experiment) {
    if (taxi_path == "UNSET") {
      std::cout << "Compression experiment requires taxi datasets." << std::endl;
    } else {
      auto taxi = enc::lib::load_taxi_table(taxi_path, taxi_scale_factor - 1);
      enc::exp::compression_experiment(taxi);
    }
  } else if (tpch_path != "UNSET") {
    std::cout << "run_number,experiment,"
                 "total_time,lt-shipdate,"
                 "geq-discount,lt-quantity,"
                 "geq-shipdate,leq-discount,agg,"
                 "verify-popcnt,verify-agg"
              << std::endl;
    auto tpch_lineitem = enc::lib::load_tpch_lineitem_table(tpch_path, db_scale_factor, false);
    enc::bench::tpch_ete_q6_byte_wrapper(tpch_lineitem, repeats);
    enc::bench::tpch_ete_q6_bit_wrapper(tpch_lineitem, repeats);
  } else if (ssb_wide_path != "UNSET") {
    std::cout << "run_number,experiment,"
                 "time,verify-popcnt"
              << std::endl;
    enc::lib::SsbWideTable ssb_wide_table = enc::lib::load_ssb_widetable(ssb_wide_path, db_scale_factor);
    enc::bench::ssb_wrapper_bit(ssb_wide_table, repeats);
    enc::bench::ssb_wrapper_byte(ssb_wide_table, repeats);
  } else if (taxi_path != "UNSET") {
    std::cout << "run_number,dataset-size,experiment,method,encoding,"
                 "time,verify-popcnt"
              << std::endl;
    auto taxi = enc::lib::load_taxi_table(taxi_path, taxi_scale_factor - 1);
    enc::bench::microbenchmarks(taxi, repeats);
  }
  return 0;
}
