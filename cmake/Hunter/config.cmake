# This file specifies additional data for the dependencies that are imported via hunter

if(APPLE)
hunter_config(Boost 
  VERSION 1.86.0
  URL "https://archives.boost.io/release/1.86.0/source/boost_1_86_0.tar.bz2"
  SHA1
  fd0d26a7d5eadf454896942124544120e3b7a38f
  CMAKE_ARGS CMAKE_CXX_FLAGS=-std=c++20 # Required for OSX, otherwise it compiles with C++98
)
else()
hunter_config(Boost 
  VERSION 1.86.0
  URL "https://archives.boost.io/release/1.86.0/source/boost_1_86_0.tar.bz2"
  SHA1
  fd0d26a7d5eadf454896942124544120e3b7a38f
)
endif()

hunter_config(Eigen VERSION 3.4.0)
hunter_config(GTest VERSION 1.11.0)
