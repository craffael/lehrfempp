cmake_minimum_required( VERSION 3.10)

include("cmake/HunterGate.cmake")

HunterGate(
    URL "https://github.com/ruslo/hunter/archive/v0.20.35.tar.gz"
    SHA1 "6e3cb4c333b76803a83c56fcbca3b0bfd9f96f27"
)

project(LehrFEMpp VERSION 1.0.0 LANGUAGES CXX)

# Check that there is no in-source build:
if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
   message(SEND_ERROR "In-source builds are not allowed.")
endif()

# Enable testing
enable_testing()

# Get boost
hunter_add_package(Boost COMPONENTS chrono timer system)
find_package(Boost CONFIG REQUIRED chrono timer system)

hunter_add_package(Boost COMPONENTS system filesystem program_options)
find_package(Boost CONFIG REQUIRED system filesystem program_options)

# Get Eigen
hunter_add_package(Eigen)
find_package(Eigen3 CONFIG REQUIRED)

# Get Google Test
hunter_add_package(GTest)
find_package(GTest CONFIG REQUIRED)

add_subdirectory(lib)
add_subdirectory(doc/doxygen)
add_subdirectory(examples)
add_subdirectory(experiments)
