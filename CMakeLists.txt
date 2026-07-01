cmake_minimum_required(VERSION 3.10)
project(Protocol)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--no-relax")

include_directories(${CMAKE_SOURCE_DIR}/mcl/include)
link_directories(${CMAKE_SOURCE_DIR}/mcl/lib)

find_library(GMP_LIB gmp)
find_library(GMPXX_LIB gmpxx)

# 编译主协议
add_executable(main zk_PIOP.cpp kzg.cpp subgroups.cpp)
target_link_libraries(main ${CMAKE_SOURCE_DIR}/mcl/lib/libmcl.a ${GMP_LIB} ${GMPXX_LIB})

# 编译 Bulletproofs 基线
add_executable(run_bulletproof bulletproof.cpp kzg.cpp subgroups.cpp)
target_link_libraries(run_bulletproof ${CMAKE_SOURCE_DIR}/mcl/lib/libmcl.a ${GMP_LIB} ${GMPXX_LIB})

# 编译 Plookup 基线
add_executable(run_plookup plookup.cpp kzg.cpp subgroups.cpp)
target_link_libraries(run_plookup ${CMAKE_SOURCE_DIR}/mcl/lib/libmcl.a ${GMP_LIB} ${GMPXX_LIB})
