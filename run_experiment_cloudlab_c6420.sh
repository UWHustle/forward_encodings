#!/bin/bash
SCRIPT=$(readlink -f "$0")
SCRIPT_PATH=$(dirname "$SCRIPT")
DB_SF=100
TAXI_SF=100
REPEATS=10
# c6420 has two sockets. bind openmp to 2nd socket.
export OMP_NUM_THREADS=32;
export GOMP_CPU_AFFINITY="1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 33 35 37 39 41 43 45 47 49 51 53 55 57 59 61 63";
echo "OMP Threads: ${OMP_NUM_THREADS}";
echo "OMP Affinity: ${GOMP_CPU_AFFINITY}";
# script assumes a user owned dir has been set up at /mydata/data
# refer to ./cloudlab_setup.sh
TPCH_PATH="/mydata/data/tpch-dbgen/build/sf${DB_SF}-lineitem.tbl";
SSB_PATH="/mydata/data/sf${DB_SF}-ssb_wide_table.tbl";
TAXI_PATH="/mydata/data/yellow_tripdata_2022-01.csv";
cd "${SCRIPT_PATH}" || exit;
python3 fetch_data.py --scale_factor="${DB_SF}" --data_path="/mydata/data/";
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
numactl --membind=1 --cpunodebind=1 -- release/encodings --repeats="${REPEATS}" --db_scale_factor="${DB_SF}" --tpch_path="${TPCH_PATH}" > "results/tpch.csv"
echo "Running benchmark ssb..."
numactl --membind=1 --cpunodebind=1 -- release/encodings --repeats="${REPEATS}" --db_scale_factor="${DB_SF}" --ssb_wide_path="${SSB_PATH}" > "results/ssb.csv"
echo "Running benchmark microbenchmark..."
numactl --membind=1 --cpunodebind=1 -- release/encodings --repeats="${REPEATS}" --taxi_scale_factor="${TAXI_SF}" --taxi_path="${TAXI_PATH}" > "results/microbench.csv"
echo "Running compression experiment..."
numactl --membind=1 --cpunodebind=1 -- release/encodings --experiment --taxi_scale_factor="${TAXI_SF}" --taxi_path="${TAXI_PATH}" > "results/experiment.csv"
