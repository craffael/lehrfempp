set(sources 
  base.h
  dereference_lambda_random_access_iterator.h
  forward_iterator.h
  forward_range.h
  invalid_type_exception.h
  lf_assert.cc
  lf_assert.h
  lf_exception.h
  random_access_iterator.h
  random_access_range.h
  ref_el.cc
  ref_el.h
  static_vars.h
  static_vars.cc
)

add_library(lf.base ${sources})
target_link_libraries(lf.base PUBLIC Eigen3::Eigen Boost::boost)
target_compile_features(lf.base PUBLIC cxx_std_17)
target_include_directories(lf.base PUBLIC
  "$<BUILD_INTERFACE:${LOCAL_INCLUDE_DIRECTORY}>"
  "$<BUILD_INTERFACE:${LOCAL_GENERATED_INCLUDE_DIRECTORY}>"
)

# Installation
set(config_install_dir "lib/cmake/lehrfempp")
set(include_install_dir "include/lf/base")

set(generated_dir "${CMAKE_CURRENT_BINARY_DIR}/generated")
set(TARGETS_EXPORT_NAME "lf.base.Targets")
set(namespace "lf::base::")

install(
    TARGETS lf.base
    EXPORT "${TARGETS_EXPORT_NAME}"
    LIBRARY DESTINATION "lib"
    ARCHIVE DESTINATION "lib"
    RUNTIME DESTINATION "bin"
    INCLUDES DESTINATION "${include_install_dir}"
)

# TODO: include only header files here
install(
    FILES ${sources}
    DESTINATION "${include_install_dir}"
)

install(
    EXPORT "${TARGETS_EXPORT_NAME}"
    NAMESPACE "${namespace}"
    DESTINATION "${config_install_dir}"
)

add_subdirectory(test)
