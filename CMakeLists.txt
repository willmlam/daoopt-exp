# CMake file for DAOOPT solver
# by Lars Otten, 2012

cmake_minimum_required(VERSION 2.6)
project(DAOOPT)

# To select the two master modes
option(DAOOPT_WORKER "Force worker binary" OFF)
option(DAOOPT_MASTER_STATIC "Force static master binary" OFF)
option(DAOOPT_MASTER_DYNAMIC "Force dynamic master binary" OFF)

# To enable static linking of the final daoopt binary
option(DAOOPT_LINK_STATIC "Link binary statically" OFF)

option(FGLP_PREPROC_ONLY "Compile FGLP preprocessor only" OFF)

# General Compiler flags
add_definitions(-Wall)

if(WIN32)
  add_definitions(-DWINDOWS)
else()
  add_definitions(-DLINUX)
endif()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS_DEBUG "-DDEBUG -O0 -g3")
set(CMAKE_C_FLAGS_DEBUG "-DDEBUG -O0 -g3")

set(CMAKE_CXX_FLAGS_RELEASE "-DDEBUG -O3")
set(CMAKE_C_FLAGS_RELEASE "-DDEBUG -O3")

# Check for C++11 support
include(CheckCXXCompilerFlag) 
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11) 
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X) 
if(COMPILER_SUPPORTS_CXX11) 
        set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -std=c++11")
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X) 
        set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -std=c++0x")
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -std=c++0x")
else() 
        message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.") 
endif() 


if(DAOOPT_LINK_STATIC)
  set(CMAKE_EXE_LINKER_FLAGS "-static")
  set(Boost_USE_STATIC_LIBS ON)
endif()

# Optional forced version: static over dynamic master, worker over both
if(DAOOPT_WORKER)
  add_definitions(-DNOTHREADS)
else()
  if(DAOOPT_MASTER_STATIC)
    add_definitions(-DPARALLEL_STATIC)
  endif()
  if(DAOOPT_MASTER_DYNAMIC AND NOT DAOOPT_MASTER_STATIC)
    add_definitions(-DPARALLEL_DYNAMIC)
    find_package(Threads)  # includes pthread
    find_package( Boost REQUIRED COMPONENTS system )
  endif()
endif()

find_package( Boost REQUIRED COMPONENTS program_options thread )
include_directories( ${Boost_INCLUDE_DIRS} )

# Included folders for header files
include_directories(
  ${DAOOPT_SOURCE_DIR}/include
  ${DAOOPT_SOURCE_DIR}/lib
)

# Libraries in the lib folder via recursive CMakeLists.txt
if (FGLP_PREPROC_ONLY)
  add_subdirectory(lib/zlib)
  add_library(Gzstream lib/gzstream.cpp)
  target_link_libraries(Gzstream Zlib)
else()
  add_subdirectory(lib)
  set(LIBS Gzstream Sls4mpe Mex Cvo)
endif()

# Main daoopt source files
set(FILES
  ./source/BestFirst.cpp
  ./source/BoundPropagator.cpp
  ./source/BoundPropagatorMaster.cpp
  ./source/BranchAndBound.cpp
  ./source/BranchAndBoundMaster.cpp
  ./source/BranchAndBoundRotate.cpp
  ./source/BranchAndBoundSampler.cpp
  ./source/CacheTable.cpp
  ./source/FGLP.cpp
  ./source/FGLPHeuristic.cpp
  ./source/FGLPMBEHybrid.cpp
  ./source/Function.cpp
  ./source/Graph.cpp
  ./source/LearningEngine.cpp
  ./source/LimitedDiscrepancy.cpp
  ./source/Main.cpp
  ./source/MBEHeuristicInstance.cpp
  ./source/MiniBucket.cpp
  ./source/MiniBucketElim.cpp
  ./source/ParallelManager.cpp
  ./source/Problem.cpp
  ./source/ProgramOptions.cpp
  ./source/Pseudotree.cpp
  ./source/Random.cpp
  ./source/PriorityFGLP.cpp
  ./source/ResidualFGLP.cpp
  ./source/Search.cpp
  ./source/SearchMaster.cpp
  ./source/SearchNode.cpp
  ./source/SigHandler.cpp
  ./source/SLSWrapper.cpp
  ./source/Statistics.cpp
  ./source/SubproblemCondor.cpp
  ./source/SubprobStats.cpp
  ./source/utils.cpp
)

set(FILES_FGLP-PREPROCESS
  ./source/CacheTable.cpp
  ./source/FGLP.cpp
  ./source/Function.cpp
  ./source/Graph.cpp
  ./source/Problem.cpp
  ./source/ProgramOptions.cpp
  ./source/Pseudotree.cpp
  ./source/Random.cpp
  ./source/utils.cpp
)

# Main executable and library dependencies
if (FGLP_PREPROC_ONLY)
  add_executable(fglpPreproc fglpPreproc.cpp ${FILES_FGLP-PREPROCESS})
  target_link_libraries(fglpPreproc Gzstream ${Boost_PROGRAM_OPTIONS_LIBRARY} ${Boost_THREAD_LIBRARY} ${Boost_SYSTEM_LIBRARY} )
else()
  add_executable(daoopt daoopt.cpp ${FILES})
  target_link_libraries(daoopt ${LIBS} ${CMAKE_THREAD_LIBS_INIT} ${Boost_PROGRAM_OPTIONS_LIBRARY} ${Boost_THREAD_LIBRARY} ${Boost_SYSTEM_LIBRARY} )
endif()
