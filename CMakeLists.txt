# CMake file for DAOOPT solver

cmake_minimum_required(VERSION 2.8)
project(daoopt_solver)

# To enable static linking of the final daoopt binary
option(DAOOPT_LINK_STATIC "Link binary statically" OFF)

if(WIN32)
  add_definitions(-DWINDOWS)
else()
  add_definitions(-DLINUX)
endif()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

if(NOT WIN32)
  set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g3 -w -fpermissive -fPIC")
  set(CMAKE_C_FLAGS_DEBUG "-O0 -g3 -fPIC")
  set(CMAKE_CXX_FLAGS_RELEASE "-O3 -w -fpermissive -fPIC")
  set(CMAKE_C_FLAGS_RELEASE "-O3 -fPIC")
else()
  set(CMAKE_CXX_FLAGS_DEBUG "-DDEBUG /Od /J")
  set(CMAKE_C_FLAGS_DEBUG "-DDEBUG /Od /J")
  set(CMAKE_CXX_FLAGS_RELEASE "/Ox /J")
  set(CMAKE_C_FLAGS_RELEASE "/Ox /J")
endif()


if(NOT WIN32)
  include(CheckCXXCompilerFlag)
  CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
  if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -std=c++11")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -std=c++11")
  else()
    message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
  endif()
endif()

if(DAOOPT_LINK_STATIC)
  set(CMAKE_EXE_LINKER_FLAGS "-static")
  set(Boost_USE_STATIC_LIBS ON)
  set(gflags_USE_STATIC_LIBS ON)
endif()

find_package(Threads)
find_package(Boost)
find_package(gflags REQUIRED)

include_directories(${gflags_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIRS})

include_directories(
  ./
  ./ARP
  ./ARP/BE
  ./ARP/CVO
  ./ARP/Problem
  ./ARP/Utils
  ./daoopt/source
  ./daoopt/lib
  ./daoopt/lib/google-utils
  ./daoopt/lib/sls4mpe
  ./mex
  ./mex/include
)

# Kalev's library, including cvo
add_subdirectory(ARP)

# daoopt
add_subdirectory(daoopt)

# Mex mplp code by Alex Ihler
add_subdirectory(mex)

# MiniSAT
add_subdirectory(miniSAT)

# Main executable
add_executable(daoopt_solver daoopt_solver.cpp
  $<TARGET_OBJECTS:ARP>
  $<TARGET_OBJECTS:daoopt>
  $<TARGET_OBJECTS:Gzstream>
  $<TARGET_OBJECTS:mex>
  $<TARGET_OBJECTS:Minisat>
  $<TARGET_OBJECTS:Sls4mpe>
  $<TARGET_OBJECTS:Zlib>)
if(DAOOPT_LINK_STATIC)
  SET_TARGET_PROPERTIES(daoopt_solver PROPERTIES LINK_SEARCH_START_STATIC 1)
  SET_TARGET_PROPERTIES(daoopt_solver PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()
target_link_libraries(daoopt_solver
  gflags
  ${Boost_THREAD_LIBRARY}
  ${Boost_SYSTEM_LIBRARY}
  ${CMAKE_THREAD_LIBS_INIT})
