#!/bin/sh -x
sudo apt-get update;
sudo apt-get --assume-yes install cmake python3-pip numactl zlib1g-dev;
python3 -m pip install Cython setuptools;
python3 -m pip install pandas duckdb pyarrow;
