cmake_minimum_required(VERSION 3.12)
project(VolatilitySmileTwoModels LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
add_executable(VolatilitySmileTwoModels
    src/main.cpp
    src/lognormal_mixture.cpp
    src/market_data.cpp
    src/multi_maturity.cpp
    src/cubic_variance.cpp)
target_include_directories(VolatilitySmileTwoModels PRIVATE src)
if(MSVC)
  target_compile_options(VolatilitySmileTwoModels PRIVATE /W4)
else()
  target_compile_options(VolatilitySmileTwoModels PRIVATE -Wall -Wextra -Wpedantic)
endif()
