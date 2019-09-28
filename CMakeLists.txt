cmake_minimum_required( VERSION 3.10)

include("cmake/functions.cmake")

HunterGate(
    URL "https://github.com/ruslo/hunter/archive/v0.23.155.tar.gz"
    SHA1 "6421d63cc4fd78b8622e863c6b1dcabce98e1216"
    LOCAL # use cmake/Hunter/config.cmake
)

project(LehrFEMpp VERSION 1.0.0 LANGUAGES CXX)

# Check that this is not an in-source build:
if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
   message(SEND_ERROR "In-source builds are not allowed.")
endif()

# CMake Options
###############################################################################
option(LF_ENABLE_TESTING "Whether the tests should be built as well" ON)
if(LF_ENABLE_TESTING) 
  # Enable testing
  enable_testing()
endif()

option(LF_BUILD_EXAMPLES "Whether the examples and experiments should be built" ON)



# Get Dependencies
# (don't forget to update cmake/Config.cmake.in !
###############################################################################
# Get boost
hunter_add_package(Boost COMPONENTS chrono timer system filesystem program_options)
find_package(Boost CONFIG REQUIRED chrono timer system filesystem program_options)

# Get Eigen
hunter_add_package(Eigen)
find_package(Eigen3 CONFIG REQUIRED)


# Get Google Test (required by test_utils modules!)
hunter_add_package(GTest)
find_package(GTest CONFIG REQUIRED)

# Add cling support
option(LF_ENABLE_CLING "if set to true, we will link with libdl so that cling works" OFF)
if("${LF_ENABLE_CLING}")
  find_library (LIBDL_PATH NAMES dl)
  set(CMAKE_CXX_STANDARD_LIBRARIES "${LIBDL_PATH}")
  if(LIBDL_PATH STREQUAL "LIBDL_PATH-NOTFOUND") 
    message(FATAL_ERROR "LF_ENABLE_CLING=On but cannot find Libdl. Consider setting LIBDL_PATH manually.")
  endif()
endif()

# Settings for installation
###############################################################################
set(LF_VERSION 0.7.3)
set(config_install_dir "lib/cmake/lehrfempp")

# Include Subdirectories
###############################################################################


add_subdirectory(lib)
add_subdirectory(doc/doxygen)
if(LF_BUILD_EXAMPLES)
  add_subdirectory(examples)
  add_subdirectory(experiments)
  add_subdirectory(projects)
endif()


# Installation
set(generated_dir "${CMAKE_CURRENT_BINARY_DIR}/generated")
set(version_config "${generated_dir}/lehrfemppConfigVersion.cmake")
set(project_config "${generated_dir}/lehrfemppConfig.cmake")

include(CMakePackageConfigHelpers)

write_basic_package_version_file("${version_config}" VERSION ${LF_VERSION} COMPATIBILITY SameMajorVersion)
configure_package_config_file(
  "cmake/Config.cmake.in"
  "${project_config}"
  INSTALL_DESTINATION "${config_install_dir}"
)

install(
    FILES "${project_config}" "${version_config}"
    DESTINATION "${config_install_dir}"
)

set(LOCAL_GENERATED_INCLUDE_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/lib)
install(EXPORT LFTargets DESTINATION ${config_install_dir} NAMESPACE LF::)



