cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

set (CMAKE_CXX_STANDARD 14)
enable_language(CXX)
enable_language(C)

set(CMAKE_VERBOSE_MAKEFILE ON)

project(Breaker)

if(DOWNLOAD_BOOST)
  project(EXTERNAL NONE)
  include (TimeEvolver/cmake/CMakeLists.txt)
  return()
else()
  project(TimeEvolver)
endif()

set (Boost_USE_STATIC_LIBS ON)

if(NOSYSTEMBOOST)
  set (Boost_NO_SYSTEM_PATHS ON)
endif()

add_subdirectory(TimeEvolver)


find_package(Boost 1.75 REQUIRED COMPONENTS program_options)


find_package(HDF5 COMPONENTS CXX)
if(HDF5_FOUND)
  add_compile_definitions(USE_HDF)
  include_directories(${HDF5_INCLUDE_DIR})
  set(_hdf5_libs hdf5 hdf5_cpp)
endif()

if(DEFINED ENV{MKLROOT})
message("MKLROOT defined")
ELSE()
include_directories(/usr/include/mkl)
ENDIF()

add_executable(Breaker quantumBreaking.cpp quantumBreakingHamiltonians.cpp)
set_target_properties(Breaker PROPERTIES OUTPUT_NAME "quantumBreaker" INSTALL_RPATH_USE_LINK_PATH TRUE)
target_include_directories(Breaker PUBLIC TimeEvolver/core TimeEvolver/helper)

target_link_libraries(Breaker PUBLIC
    ${BLAS_LIBRARIES}
    ${Boost_LIBRARIES}
    dl pthread
    Helper TimeEvolver) 

if(HDF5_FOUND)
    target_link_libraries(Breaker PUBLIC
    ${HDF5_LIBRARIES}) 
endif()