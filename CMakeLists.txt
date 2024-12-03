<<<<<<< HEAD
cmake_minimum_required(VERSION 3.10)
project(Protocol)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--no-relax")

include_directories(${CMAKE_SOURCE_DIR}/mcl/include)
link_directories(${CMAKE_SOURCE_DIR}/mcl/lib)

find_library(GMP_LIB gmp)
find_library(GMPXX_LIB gmpxx)

add_executable(main zk_PIOP.cpp kzg.cpp subgroups.cpp)

target_link_libraries(main mcl mclbn256 ${GMP_LIB} ${GMPXX_LIB})
=======
cmake_minimum_required(VERSION 3.10)
project(Protocol)

set(CMAKE_CXX_STANDARD 17)

include_directories(${CMAKE_SOURCE_DIR}/mcl/include)
link_directories(${CMAKE_SOURCE_DIR}/mcl/lib)

find_library(GMP_LIB gmp)
find_library(GMPXX_LIB gmpxx)

add_executable(main zk_PIOP.cpp kzg.cpp subgroups.cpp)

target_link_libraries(main mcl ${GMP_LIB} ${GMPXX_LIB})
>>>>>>> ea97283b7d0d33a96163f269cae4f9e92a8490e0
