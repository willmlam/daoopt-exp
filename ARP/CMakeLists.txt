include_directories(../minisat)

add_library(ARP OBJECT
  ./BE/MiniBucket.cpp
  ./BE/Bucket.cpp
  ./BE/MBEworkspace.cpp
  ./CVO/Graph.cpp
  ./CVO/Graph_MinFillOrderComputation.cpp
  ./CVO/VariableOrderComputation.cpp
  ./CVO/Graph_RemoveRedundantFillEdges.cpp
  ./Problem/Problem.cpp
  ./Problem/Globals.cpp
  ./Problem/Workspace.cpp
  ./Problem/Function.cpp
  ./Utils/AVLtreeSimple.cpp
  ./Utils/Mutex.cpp
  ./Utils/MiscUtils.cpp
  ./Utils/FnExecutionThread.cpp
  ./Utils/Sort.cxx
)
