### Fast Predicate-based Columnar Scans Using Forward Encodings

Submitted to PVLDB Vol. 16.

### Notice

This codebase uses AVX2 instructions, and may not work when run using a machine
with a non-Intel CPU.

### Instructions

Experiments were performed using a CloudLab c6420 machine, running
a standard Ubuntu 22 LTS distribution. 
This machine was configured to have 100GB of storage mounted to `/mydata/`.

Dependencies were automatically installed using the `clooudlab_setup.sh` 
script.
Experiments were run by executing the `run_experiment_cloudlab_c6420.sh` script.
The executor script wraps around both `fetch_data.py` and the compiled experiment
binary.

These scripts have been adapted to `simple_setup.sh` 
and `simple_run_experiment_sf1.sh`,
which can be run on any Intel-based Ubuntu 22 machine.
The script assumes GCC-11 is located at `/usr/bin/gcc-11`.
Data will be fetched to a local `data/` directory, experiment binary
will be built in `release/`, and results will be placed in `results/`.
The experiment is run at scale factor 1.

To run the experiments:
```
$ ./simple_setup.sh
$ ./simple_run_experiment_sf1.sh
```

Experiment results are in the form of 4 CSV files:

`experiment.csv`: The compression experiment.

`microbench.csv`: The microbenchmarks using the Taxi dataset.

`ssb.csv`: The SSB results.

`tpch.csv`: The TPC-H Q6 results.

All times reported in the CSV files are in nanoseconds.
