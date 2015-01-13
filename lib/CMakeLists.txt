# Zlib library
add_subdirectory(zlib)

# SLS4MPE code by Frank Hutter
add_subdirectory(sls4mpe)

# GzStream library (needs gzip)
add_library(Gzstream gzstream.cpp)
target_link_libraries(Gzstream Zlib)

# Mex mplp code by Alex Ihler
add_subdirectory(mex)

# Kalev's library, including cvo
add_subdirectory(ARP)

# MiniSAT
add_subdirectory(miniSAT)

# google3-like utilities
#add_subdirectory(google-utils)
