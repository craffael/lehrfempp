add_subdirectory(interface)
add_subdirectory(geom)
add_subdirectory(test_utils)

if(OpenCASCADE_FOUND)
  add_subdirectory(occt)
endif()