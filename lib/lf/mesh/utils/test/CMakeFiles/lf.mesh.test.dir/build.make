# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.11

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /hg/u/magina/Documents/lehrfempp/lib/lf/mesh

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /hg/u/magina/Documents/lehrfempp/lib/lf/mesh/utils

# Include any dependencies generated for this target.
include test/CMakeFiles/lf.mesh.test.dir/depend.make

# Include the progress variables for this target.
include test/CMakeFiles/lf.mesh.test.dir/progress.make

# Include the compile flags for this target's objects.
include test/CMakeFiles/lf.mesh.test.dir/flags.make

test/CMakeFiles/lf.mesh.test.dir/tp_triag_mesh_builder_tests.o: test/CMakeFiles/lf.mesh.test.dir/flags.make
test/CMakeFiles/lf.mesh.test.dir/tp_triag_mesh_builder_tests.o: ../test/tp_triag_mesh_builder_tests.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/hg/u/magina/Documents/lehrfempp/lib/lf/mesh/utils/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object test/CMakeFiles/lf.mesh.test.dir/tp_triag_mesh_builder_tests.o"
	cd /hg/u/magina/Documents/lehrfempp/lib/lf/mesh/utils/test && /usr/lib64/ccache/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/lf.mesh.test.dir/tp_triag_mesh_builder_tests.o -c /hg/u/magina/Documents/lehrfempp/lib/lf/mesh/test/tp_triag_mesh_builder_tests.cc

test/CMakeFiles/lf.mesh.test.dir/tp_triag_mesh_builder_tests.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/lf.mesh.test.dir/tp_triag_mesh_builder_tests.i"
	cd /hg/u/magina/Documents/lehrfempp/lib/lf/mesh/utils/test && /usr/lib64/ccache/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /hg/u/magina/Documents/lehrfempp/lib/lf/mesh/test/tp_triag_mesh_builder_tests.cc > CMakeFiles/lf.mesh.test.dir/tp_triag_mesh_builder_tests.i

test/CMakeFiles/lf.mesh.test.dir/tp_triag_mesh_builder_tests.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/lf.mesh.test.dir/tp_triag_mesh_builder_tests.s"
	cd /hg/u/magina/Documents/lehrfempp/lib/lf/mesh/utils/test && /usr/lib64/ccache/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /hg/u/magina/Documents/lehrfempp/lib/lf/mesh/test/tp_triag_mesh_builder_tests.cc -o CMakeFiles/lf.mesh.test.dir/tp_triag_mesh_builder_tests.s

# Object files for target lf.mesh.test
lf_mesh_test_OBJECTS = \
"CMakeFiles/lf.mesh.test.dir/tp_triag_mesh_builder_tests.o"

# External object files for target lf.mesh.test
lf_mesh_test_EXTERNAL_OBJECTS =

test/lf.mesh.test: test/CMakeFiles/lf.mesh.test.dir/tp_triag_mesh_builder_tests.o
test/lf.mesh.test: test/CMakeFiles/lf.mesh.test.dir/build.make
test/lf.mesh.test: hybrid2d/liblf.mesh.hybrid2d.a
test/lf.mesh.test: hybrid2dp/liblf.mesh.hybrid2dp.a
test/lf.mesh.test: utils/liblf.mesh.utils.a
test/lf.mesh.test: test_utils/liblf.mesh.test_utils.a
test/lf.mesh.test: liblf.mesh.a
test/lf.mesh.test: test/CMakeFiles/lf.mesh.test.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/hg/u/magina/Documents/lehrfempp/lib/lf/mesh/utils/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable lf.mesh.test"
	cd /hg/u/magina/Documents/lehrfempp/lib/lf/mesh/utils/test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/lf.mesh.test.dir/link.txt --verbose=$(VERBOSE)
	cd /hg/u/magina/Documents/lehrfempp/lib/lf/mesh/utils/test && /usr/bin/cmake -D TEST_TARGET=lf.mesh.test -D TEST_EXECUTABLE=/hg/u/magina/Documents/lehrfempp/lib/lf/mesh/utils/test/lf.mesh.test -D TEST_EXECUTOR= -D TEST_WORKING_DIR=/hg/u/magina/Documents/lehrfempp/lib/lf/mesh/utils/test -D TEST_EXTRA_ARGS= -D TEST_PROPERTIES= -D TEST_PREFIX= -D TEST_SUFFIX= -D NO_PRETTY_TYPES=FALSE -D NO_PRETTY_VALUES=FALSE -D TEST_LIST=lf.mesh.test_TESTS -D CTEST_FILE=/hg/u/magina/Documents/lehrfempp/lib/lf/mesh/utils/test/lf.mesh.test[1]_tests.cmake -D TEST_DISCOVERY_TIMEOUT=5 -P /usr/share/cmake/Modules/GoogleTestAddTests.cmake

# Rule to build all files generated by this target.
test/CMakeFiles/lf.mesh.test.dir/build: test/lf.mesh.test

.PHONY : test/CMakeFiles/lf.mesh.test.dir/build

test/CMakeFiles/lf.mesh.test.dir/clean:
	cd /hg/u/magina/Documents/lehrfempp/lib/lf/mesh/utils/test && $(CMAKE_COMMAND) -P CMakeFiles/lf.mesh.test.dir/cmake_clean.cmake
.PHONY : test/CMakeFiles/lf.mesh.test.dir/clean

test/CMakeFiles/lf.mesh.test.dir/depend:
	cd /hg/u/magina/Documents/lehrfempp/lib/lf/mesh/utils && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /hg/u/magina/Documents/lehrfempp/lib/lf/mesh /hg/u/magina/Documents/lehrfempp/lib/lf/mesh/test /hg/u/magina/Documents/lehrfempp/lib/lf/mesh/utils /hg/u/magina/Documents/lehrfempp/lib/lf/mesh/utils/test /hg/u/magina/Documents/lehrfempp/lib/lf/mesh/utils/test/CMakeFiles/lf.mesh.test.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : test/CMakeFiles/lf.mesh.test.dir/depend
