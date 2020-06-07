cmake_minimum_required(VERSION 3.10)

### Hunter package manager ###
include("cmake/HunterGate.cmake")
HunterGate(
    URL "https://github.com/cpp-pm/hunter/archive/v0.23.249.tar.gz"
    SHA1 "d45d77d8bba9da13e9290a180e0477e90accd89b"
	LOCAL # use cmake/Hunter/config.cmake
)

project(fisherkpp)

message("Project source dir =" "${PROJECT_SOURCE_DIR}")
message("Project binary dir =" "${PROJECT_BINARY_DIR}")

# Check that this is not an in-source build:
if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
   message(SEND_ERROR "In-source builds are not allowed.")
endif()

# Setup compiler
set(CMAKE_CXX_STANDARD 17)

# Enable warnings
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG(-Wall COMPILER_SUPPORTS_WALL)
CHECK_CXX_COMPILER_FLAG(-Wextra COMPILER_SUPPORTS_WEXTRA)
CHECK_CXX_COMPILER_FLAG(-Wconversion COMPILER_SUPPORTS_WCONVERSION)

add_library(warning_flags INTERFACE)
if(COMPILER_SUPPORTS_WALL)
  target_compile_options(warning_flags INTERFACE -Wall)
endif()
if(COMPILER_SUPPORTS_WEXTRA)
  target_compile_options(warning_flags INTERFACE -Wextra)
endif()
if(COMPILER_SUPPORTS_WCONVERSION)
  target_compile_options(warning_flags INTERFACE -Wconversion)
endif()

# Download dependencies
###############################################################################
# Get boost
hunter_add_package(Boost COMPONENTS chrono timer system filesystem program_options)
find_package(Boost CONFIG REQUIRED chrono timer system filesystem program_options)

# Get Eigen
hunter_add_package(Eigen)
find_package(Eigen3 CONFIG REQUIRED)

# Get LehrFEMpp
hunter_add_package(lehrfempp)
find_package(lehrfempp CONFIG REQUIRED)
###############################################################################

# Set sources
set(SOURCES 
	${DIR}/strangsplitting.h
	${DIR}/strangsplitting.cc
	${DIR}/strangsplitting_main.cc
	${DIR}/modelproblem_island_main.cc
	${DIR}/modelproblem_circle_main.cc
	${DIR}/modelproblem_threecircles_main.cc
	${DIR}/convergence_meshwidth_main.cc
	${DIR}/convergence_tau_main.cc
	${DIR}/convergence_linked_main.cc
	)

# Define executable
add_executable(splitting strangsplitting_main.cc)
add_executable(testIsland modelproblem_island_main.cc)
add_executable(testCircle modelproblem_circle_main.cc)
add_executable(test3Circles modelproblem_threecircles_main.cc)
add_executable(convergenceH convergence_meshwidth_main.cc)
add_executable(convergenceTau convergence_tau_main.cc)
add_executable(convergenceLinked convergence_linked_main.cc)

# link executable with required modules

target_link_libraries(splitting
					Boost::filesystem
					Boost::program_options
					Boost::system
					Eigen3::Eigen
					LF::lf.assemble
					LF::lf.base
					LF::lf.geometry
					LF::lf.io
					LF::lf.mesh
					LF::lf.mesh.hybrid2d
					LF::lf.mesh.test_utils
					LF::lf.mesh.utils
					LF::lf.refinement
					LF::lf.uscalfe
					)

target_link_libraries(testIsland
					Boost::filesystem
					Boost::program_options
					Boost::system
					Eigen3::Eigen
					LF::lf.assemble
					LF::lf.base
					LF::lf.geometry
					LF::lf.io
					LF::lf.mesh
					LF::lf.mesh.hybrid2d
					LF::lf.mesh.test_utils
					LF::lf.mesh.utils
					LF::lf.refinement
					LF::lf.uscalfe
					)

target_link_libraries(testCircle
					Boost::filesystem
					Boost::program_options
					Boost::system
					Eigen3::Eigen
					LF::lf.assemble
					LF::lf.base
					LF::lf.geometry
					LF::lf.io
					LF::lf.mesh
					LF::lf.mesh.hybrid2d
					LF::lf.mesh.test_utils
					LF::lf.mesh.utils
					LF::lf.refinement
					LF::lf.uscalfe
					)

target_link_libraries(test3Circles
					Boost::filesystem
					Boost::program_options
					Boost::system
					Eigen3::Eigen
					LF::lf.assemble
					LF::lf.base
					LF::lf.geometry
					LF::lf.io
					LF::lf.mesh
					LF::lf.mesh.hybrid2d
					LF::lf.mesh.test_utils
					LF::lf.mesh.utils
					LF::lf.refinement
					LF::lf.uscalfe
					)

target_link_libraries(convergenceH
					Boost::filesystem
					Boost::program_options
					Boost::system
					Eigen3::Eigen
					LF::lf.assemble
					LF::lf.base
					LF::lf.geometry
					LF::lf.io
					LF::lf.mesh
					LF::lf.mesh.hybrid2d
					LF::lf.mesh.test_utils
					LF::lf.mesh.utils
					LF::lf.refinement
					LF::lf.uscalfe
					)

target_link_libraries(convergenceTau
					Boost::filesystem
					Boost::program_options
					Boost::system
					Eigen3::Eigen
					LF::lf.assemble
					LF::lf.base
					LF::lf.geometry
					LF::lf.io
					LF::lf.mesh
					LF::lf.mesh.hybrid2d
					LF::lf.mesh.test_utils
					LF::lf.mesh.utils
					LF::lf.refinement
					LF::lf.uscalfe
					)

target_link_libraries(convergenceLinked
					Boost::filesystem
					Boost::program_options
					Boost::system
					Eigen3::Eigen
					LF::lf.assemble
					LF::lf.base
					LF::lf.geometry
					LF::lf.io
					LF::lf.mesh
					LF::lf.mesh.hybrid2d
					LF::lf.mesh.test_utils
					LF::lf.mesh.utils
					LF::lf.refinement
					LF::lf.uscalfe
					)

add_custom_target(splitting_run COMMAND splitting)
add_custom_target(testIsland_run COMMAND testIsland)
add_custom_target(testCircle_run COMMAND testCircle)
add_custom_target(test3Circles_run COMMAND test3Circles)
add_custom_target(convergenceH_run COMMAND convergenceH)
add_custom_target(convergenceTau_run COMMAND convergenceTau)
add_custom_target(convergenceLinked_run COMMAND convergenceLinked)
