cmake_minimum_required(VERSION 3.10)
project(Protocol)

set(CMAKE_CXX_STANDARD 17)

include_directories(${CMAKE_SOURCE_DIR}/mcl/include)
link_directories(${CMAKE_SOURCE_DIR}/mcl/lib)

find_library(GMP_LIB gmp)
find_library(GMPXX_LIB gmpxx)

add_executable(main zk_PIOP.cpp kzg.cpp subgroups.cpp)

target_link_libraries(main mcl ${GMP_LIB} ${GMPXX_LIB})
