import argparse
import subprocess
import urllib.request
import os

import duckdb
import pandas


def dbgen(args):
    subprocess.run(
        f'''
        mkdir -p {args.data_path}
        ''',
        shell=True
    )
    fpfx = f'sf{args.scale_factor}-'
    print('Fetching data...')
    if not args.skip_tpch:
        if not os.path.exists(args.data_path + "/tpch-dbgen"):
            print('Fetching and building TPC-H data generator.')
            subprocess.run(
                f'''
                cd {args.data_path}
                git clone https://github.com/eyalroz/tpch-dbgen.git
                ''',
                shell=True)
            with open(f'{args.data_path}/tpch-dbgen/CMakeLists.txt', 'r') as file:
                filedata = file.read()
                filedata = filedata.replace("set_property(TARGET dbgen qgen APPEND PROPERTY C_STANDARD 99)",
                                            "set_target_properties(dbgen qgen PROPERTIES C_STANDARD 99)")
                filedata = filedata.replace("set_property(TARGET dbgen qgen APPEND PROPERTY C_EXTENSIONS OFF)",
                                            "set_target_properties(dbgen qgen PROPERTIES C_EXTENSIONS OFF)")
            with open(f'{args.data_path}/tpch-dbgen/CMakeLists.txt', 'w') as file:
                file.write(filedata)
            subprocess.run(
                f'''
                cd {args.data_path}
                mkdir -p tpch-dbgen/build
                cmake -Btpch-dbgen/build -Htpch-dbgen -DCMAKE_BUILD_TYPE=Release
                cmake --build tpch-dbgen/build
                ''',
                shell=True)
        else:
            print(f'TPC-H data gen found, skipping...')
        if not os.path.exists(f'{args.data_path}/tpch-dbgen/build/{fpfx}lineitem.tbl'):
            print(f'Generating TPCH data for scale factor {args.scale_factor}.')
            subprocess.run(
                f'''
                cd {args.data_path}/tpch-dbgen/build
                rm -f customer.tbl date.tbl lineorder.tbl part.tbl supplier.tbl
                ./dbgen -b ../dists.dss -s {args.scale_factor}
                for filename in *.tbl; do mv "$filename" "{fpfx}${{filename}}"; done;
                ''',
                shell=True)
        else:
            print(f'{fpfx}TPC-H found, skipping...')
    if not args.skip_ssb:
        if not os.path.exists(args.data_path + "/ssb-dbgen"):
            print('Fetching and building SSB data generator.')
            subprocess.run(
                f'''
                cd {args.data_path}
                git clone https://github.com/eyalroz/ssb-dbgen.git
                mkdir -p ssb-dbgen/build
                cmake -Bssb-dbgen/build -Hssb-dbgen -DCMAKE_BUILD_TYPE=Release
                cmake --build ssb-dbgen/build
                ''',
                shell=True)
        else:
            print(f'SSB data gen found, skipping...')
        if not os.path.exists(f'{args.data_path}/ssb-dbgen/build/{fpfx}lineorder.tbl'):
            print(f'Generating SSB data for scale factor {args.scale_factor}.')
            subprocess.run(
                f'''
                cd {args.data_path}/ssb-dbgen/build
                rm -f customer.tbl lineitem.tbl nation.tbl orders.tbl part.tbl partsupp.tbl region.tbl supplier.tbl
                ./dbgen -b ../dists.dss -s {args.scale_factor}
                for filename in *.tbl; do mv "$filename" "{fpfx}${{filename}}"; done;
                ''',
                shell=True)
        else:
            print(f'{fpfx}SSB found, skipping...')
        if args.wide_ssb:
            if not os.path.exists(f'{args.data_path}/{fpfx}ssb_wide_table.tbl'):
                print('Generating SSB WideTable with DuckDB.')
                cursor = duckdb.connect()
                with open('ssb_wide_table.sql', 'r') as file:
                    filedata = file.read()
                    filedata = filedata.replace("?1", f'{args.data_path}/')
                    filedata = filedata.replace("?2", f'{fpfx}')
                    cursor.execute(filedata)
            else:
                print(f'{fpfx}SSB Widetable found, skipping...')
    if not args.skip_taxi:
        taxi_link = "https://d37ci6vzurychx.cloudfront.net/trip-data/yellow_tripdata_2022-01.parquet"
        taxi_file = args.data_path + "/yellow_tripdata_2022-01.parquet"
        taxi_csv = args.data_path + "/yellow_tripdata_2022-01.csv"
        if not os.path.exists(taxi_csv):
            if not os.path.exists(taxi_file):
                print("Downloading taxi dataset...")
                urllib.request.urlretrieve(taxi_link, taxi_file)
            print("Converting parquet to csv...")
            taxi_df = pandas.read_parquet(taxi_file)
            taxi_df.to_csv(taxi_csv, index=False)
        else:
            print(f'Taxi data found, skipping...')
    print('Done')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate SSB data.')
    parser.add_argument('--scale_factor', type=int, default=1)
    parser.add_argument('--data_path', type=str, default="res")
    parser.add_argument('--skip_tpch', dest='skip_tpch', action='store_true', default=False)
    parser.add_argument('--skip_ssb', dest='skip_ssb', action='store_true', default=False)
    parser.add_argument('--wide_ssb', dest='wide_ssb', action='store_true', default=True)
    parser.add_argument('--skip_taxi', dest='skip_taxi', action='store_true', default=False)
    parsed_args = parser.parse_args()
    dbgen(parsed_args)
