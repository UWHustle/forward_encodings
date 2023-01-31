cmake_minimum_required(VERSION 3.17)
project(encodings)

set(EVENT_DISABLE_TESTS ON)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -fopenmp -Wno-error -lpthread -std=c++20")
set(CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} -march=native -fopenmp -Wno-error -lpthread -std=c++20")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -DNDEBUG")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -DNDEBUG")

find_package(Threads REQUIRED)
find_package(OpenMP REQUIRED)
find_package(ZLIB)

set(CMAKE_BUILD_TYPE_TEMP ${CMAKE_BUILD_TYPE})
set(CMAKE_BUILD_TYPE Release)
include(FetchContent)
FetchContent_Declare(
        csv
        GIT_REPOSITORY https://github.com/vincentlaucsb/csv-parser
        GIT_TAG 2.1.3
)
FetchContent_MakeAvailable(csv)
FetchContent_Declare(
        cxxopts
        GIT_REPOSITORY https://github.com/jarro2783/cxxopts
        GIT_TAG v3.0.0
)
FetchContent_MakeAvailable(cxxopts)
FetchContent_Declare(
        date
        GIT_REPOSITORY https://github.com/HowardHinnant/date
        GIT_TAG v3.0.1
)
FetchContent_MakeAvailable(date)
set(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE_TEMP})

target_include_directories(csv INTERFACE
        ${csv_SOURCE_DIR}/include/)

add_executable(encodings
        main.cpp
        encodings.h
        stratified_storage.h
        benchmarks/bench_util.h
        benchmarks/tpch-ete.h
        benchmarks/microbenchmarks.h
        benchmarks/ssb-scans_bit.h
        benchmarks/ssb-scans_byte.h
        experiments/compression.h
        wrappers/wrapper_utilities.h wrappers/wrapper_utilities.cpp
        wrappers/ssb_wrapper.h wrappers/ssb_wrapper.cpp
        wrappers/taxi_wrapper.h wrappers/taxi_wrapper.cpp
        wrappers/tpch_wrapper.h wrappers/tpch_wrapper.cpp
        )
target_link_libraries(encodings
        Threads::Threads
        OpenMP::OpenMP_CXX
        cxxopts
        csv
        date
        ZLIB::ZLIB
        )
