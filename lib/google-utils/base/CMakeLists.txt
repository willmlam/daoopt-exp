include_directories(
  ../
)

add_library(google_base
    callback.cc
    hash.cc
    logging.cc
    notification.cc
    stringprintf.cc
    strtoint.cc
    sysinfo.cc
    timer.cc
    walltime.cc
)
