include_directories(./include)

add_library(mex OBJECT
./src/Factor.cpp
./src/graphmodel.cpp
./src/factorgraph.cpp
./mex-workspace.cpp
)
