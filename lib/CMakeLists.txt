# Zlib library
add_subdirectory(zlib)

# SLS4MPE code by Frank Hutter
add_subdirectory(sls4mpe)

# GzStream library (needs gzip)
add_library(Gzstream gzstream.cpp)
target_link_libraries(Gzstream Zlib)

# Mex mplp code by Alex Ihler
add_subdirectory(mex)

# Kalev's ordering code
add_subdirectory(cvo)
