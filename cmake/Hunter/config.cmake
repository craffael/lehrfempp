# This file specifies additional data for the dependencies that are imported via hunter

if(APPLE)
hunter_config(Boost 
  VERSION 1.86.0
  CMAKE_ARGS CMAKE_CXX_FLAGS=-std=c++20 # Required for OSX, otherwise it compiles with C++98
)
else()
hunter_config(Boost 
  VERSION 1.86.0
)
endif()

hunter_config(Eigen VERSION 3.4.0)
hunter_config(GTest VERSION 1.15.2)
