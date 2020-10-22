include(git_info.cmake)

set(sources
  base.h
  comm.h
  comm.cc
  eigen_tools.h
  invalid_type_exception.h
  lf_assert.cc
  lf_assert.h
  lf_exception.h
  predicate_true.h
  ref_el.cc
  ref_el.h
  span.h
)

lf_add_library(lf.base ${sources})
target_link_libraries(lf.base PUBLIC Eigen3::Eigen Boost::boost Boost::program_options)

if(MSVC)
  if(${MSVC_VERSION} GREATER_EQUAL 1915) 
    # You must acknowledge that you understand MSVC resolved a byte alignment issue in this compiler
    # We get this due to using Eigen objects and allocating those objects with make_shared
    target_compile_definitions( lf.base PUBLIC _ENABLE_EXTENDED_ALIGNED_STORAGE )
  endif()
endif()

if(LF_ENABLE_TESTING)
  add_subdirectory(test)
endif()


