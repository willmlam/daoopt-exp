find_package(Threads)
find_package( Boost REQUIRED COMPONENTS thread )
include_directories(${Boost_INCLUDE_DIRS})

# Included folders for header files
include_directories(
  ./source
  ./lib
  ./lib/google-utils
)


add_subdirectory(lib)

# Main daoopt source files
add_library(daoopt OBJECT
  #  ./source/BestFirst.cpp
  ./source/BoundPropagator.cpp
  ./source/BoundPropagatorMaster.cpp
  ./source/BranchAndBound.cpp
  ./source/BranchAndBoundMaster.cpp
  ./source/BranchAndBoundRotate.cpp
  ./source/BranchAndBoundSampler.cpp
  ./source/CacheTable.cpp
  ./source/DaooptInterface.cpp
  ./source/FGLP.cpp
  ./source/FGLPHeuristic.cpp
  ./source/FGLPMBEHybrid.cpp
  ./source/Function.cpp
  ./source/Graph.cpp
  ./source/LearningEngine.cpp
  ./source/LimitedDiscrepancy.cpp
  ./source/Main.cpp
  ./source/MiniBucket.cpp
  ./source/MiniBucketElim.cpp
  ./source/MiniBucketElimInc.cpp
  ./source/MiniBucketElimLH.cpp
  ./source/MiniBucketElimLHsubtree.cpp
  ./source/ParallelManager.cpp
  ./source/PriorityFGLP.cpp
  ./source/Problem.cpp
  ./source/ProgramOptions.cpp
  ./source/Pseudotree.cpp
  ./source/Random.cpp
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
