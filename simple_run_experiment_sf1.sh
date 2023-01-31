#!/bin/bash
SCRIPT=$(readlink -f "$0")
SCRIPT_PATH=$(dirname "$SCRIPT")
DB_SF=1
TAXI_SF=1
REPEATS=10
TPCH_PATH="data/tpch-dbgen/build/sf${DB_SF}-lineitem.tbl";
SSB_PATH="data/sf${DB_SF}-ssb_wide_table.tbl";
TAXI_PATH="data/yellow_tripdata_2022-01.csv";
cd "${SCRIPT_PATH}" || exit;
python3 fetch_data.py --scale_factor="${DB_SF}" --data_path="data/";
if [ ! -d "release" ]; then
  mkdir release;
  cd release || exit;
  cmake -DCMAKE_C_COMPILER=/usr/bin/gcc-11 -DCMAKE_CXX_COMPILER=/usr/bin/g++-11 -DCMAKE_BUILD_TYPE=Release ../src;
  make CC=/usr/bin/gcc-11 CXX=/usr/local/bin/g++-11 -j"$(nproc)";
else
  cd release || exit;
fi;
cd .. || exit;
if [ ! -d "results" ]; then
  mkdir results;
fi;
# c6420 has two sockets. numactl to second socket
echo "Running benchmark tpch...";
release/encodings --repeats="${REPEATS}" --db_scale_factor="${DB_SF}" --tpch_path="${TPCH_PATH}" > "results/tpch.csv"
echo "Running benchmark ssb..."
release/encodings --repeats="${REPEATS}" --db_scale_factor="${DB_SF}" --ssb_wide_path="${SSB_PATH}" > "results/ssb.csv"
echo "Running benchmark microbenchmark..."
release/encodings --repeats="${REPEATS}" --taxi_scale_factor="${TAXI_SF}" --taxi_path="${TAXI_PATH}" > "results/microbench.csv"
echo "Running compression experiment..."
release/encodings --experiment --taxi_scale_factor="${TAXI_SF}" --taxi_path="${TAXI_PATH}" > "results/experiment.csv"
