# Setup filtering of doxygen source files so that we can turn Latex references 
# from NPDELrefs.aux into section numbers in the doxygen doc.

# C++ script which is executed for every doxygen source file:
set(CMAKE_CXX_FLAGS_DEBUG ${CMAKE_CXX_FLAGS_RELEASE}) # compile in release mode, even if building debug (for efficiency)
set(CMAKE_CXX_FLAGS_MINSIZEREL ${CMAKE_CXX_FLAGS_RELEASE})

set(SOURCES
  filters.cpp
  gh_edit_filter.cpp
  link_filter.cpp
)

add_executable(doxygen_filter ${SOURCES})

target_link_libraries(doxygen_filter PUBLIC Boost::boost)
target_link_libraries(doxygen_filter PUBLIC Boost::program_options)

set_target_properties(doxygen_filter PROPERTIES 
  FOLDER "doc"
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}$<0:hello>"
  RUNTIME_OUTPUT_NAME "doxygen_filter"
  CXX_STANDARD 17
)
