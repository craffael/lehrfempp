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
include utils/CMakeFiles/lf.mesh.utils.dir/depend.make

# Include the progress variables for this target.
include utils/CMakeFiles/lf.mesh.utils.dir/progress.make

# Include the compile flags for this target's objects.
include utils/CMakeFiles/lf.mesh.utils.dir/flags.make

utils/CMakeFiles/lf.mesh.utils.dir/print_info.o: utils/CMakeFiles/lf.mesh.utils.dir/flags.make
utils/CMakeFiles/lf.mesh.utils.dir/print_info.o: print_info.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/hg/u/magina/Documents/lehrfempp/lib/lf/mesh/utils/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object utils/CMakeFiles/lf.mesh.utils.dir/print_info.o"
	cd /hg/u/magina/Documents/lehrfempp/lib/lf/mesh/utils/utils && /usr/lib64/ccache/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/lf.mesh.utils.dir/print_info.o -c /hg/u/magina/Documents/lehrfempp/lib/lf/mesh/utils/print_info.cc

utils/CMakeFiles/lf.mesh.utils.dir/print_info.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/lf.mesh.utils.dir/print_info.i"
	cd /hg/u/magina/Documents/lehrfempp/lib/lf/mesh/utils/utils && /usr/lib64/ccache/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /hg/u/magina/Documents/lehrfempp/lib/lf/mesh/utils/print_info.cc > CMakeFiles/lf.mesh.utils.dir/print_info.i

utils/CMakeFiles/lf.mesh.utils.dir/print_info.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/lf.mesh.utils.dir/print_info.s"
	cd /hg/u/magina/Documents/lehrfempp/lib/lf/mesh/utils/utils && /usr/lib64/ccache/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /hg/u/magina/Documents/lehrfempp/lib/lf/mesh/utils/print_info.cc -o CMakeFiles/lf.mesh.utils.dir/print_info.s

utils/CMakeFiles/lf.mesh.utils.dir/write_matlab.o: utils/CMakeFiles/lf.mesh.utils.dir/flags.make
utils/CMakeFiles/lf.mesh.utils.dir/write_matlab.o: write_matlab.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/hg/u/magina/Documents/lehrfempp/lib/lf/mesh/utils/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object utils/CMakeFiles/lf.mesh.utils.dir/write_matlab.o"
	cd /hg/u/magina/Documents/lehrfempp/lib/lf/mesh/utils/utils && /usr/lib64/ccache/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/lf.mesh.utils.dir/write_matlab.o -c /hg/u/magina/Documents/lehrfempp/lib/lf/mesh/utils/write_matlab.cc

utils/CMakeFiles/lf.mesh.utils.dir/write_matlab.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/lf.mesh.utils.dir/write_matlab.i"
	cd /hg/u/magina/Documents/lehrfempp/lib/lf/mesh/utils/utils && /usr/lib64/ccache/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /hg/u/magina/Documents/lehrfempp/lib/lf/mesh/utils/write_matlab.cc > CMakeFiles/lf.mesh.utils.dir/write_matlab.i

utils/CMakeFiles/lf.mesh.utils.dir/write_matlab.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/lf.mesh.utils.dir/write_matlab.s"
	cd /hg/u/magina/Documents/lehrfempp/lib/lf/mesh/utils/utils && /usr/lib64/ccache/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /hg/u/magina/Documents/lehrfempp/lib/lf/mesh/utils/write_matlab.cc -o CMakeFiles/lf.mesh.utils.dir/write_matlab.s

utils/CMakeFiles/lf.mesh.utils.dir/write_tikz.o: utils/CMakeFiles/lf.mesh.utils.dir/flags.make
utils/CMakeFiles/lf.mesh.utils.dir/write_tikz.o: write_tikz.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/hg/u/magina/Documents/lehrfempp/lib/lf/mesh/utils/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object utils/CMakeFiles/lf.mesh.utils.dir/write_tikz.o"
	cd /hg/u/magina/Documents/lehrfempp/lib/lf/mesh/utils/utils && /usr/lib64/ccache/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/lf.mesh.utils.dir/write_tikz.o -c /hg/u/magina/Documents/lehrfempp/lib/lf/mesh/utils/write_tikz.cc

utils/CMakeFiles/lf.mesh.utils.dir/write_tikz.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/lf.mesh.utils.dir/write_tikz.i"
	cd /hg/u/magina/Documents/lehrfempp/lib/lf/mesh/utils/utils && /usr/lib64/ccache/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /hg/u/magina/Documents/lehrfempp/lib/lf/mesh/utils/write_tikz.cc > CMakeFiles/lf.mesh.utils.dir/write_tikz.i

utils/CMakeFiles/lf.mesh.utils.dir/write_tikz.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/lf.mesh.utils.dir/write_tikz.s"
	cd /hg/u/magina/Documents/lehrfempp/lib/lf/mesh/utils/utils && /usr/lib64/ccache/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /hg/u/magina/Documents/lehrfempp/lib/lf/mesh/utils/write_tikz.cc -o CMakeFiles/lf.mesh.utils.dir/write_tikz.s

# Object files for target lf.mesh.utils
lf_mesh_utils_OBJECTS = \
"CMakeFiles/lf.mesh.utils.dir/print_info.o" \
"CMakeFiles/lf.mesh.utils.dir/write_matlab.o" \
"CMakeFiles/lf.mesh.utils.dir/write_tikz.o"

# External object files for target lf.mesh.utils
lf_mesh_utils_EXTERNAL_OBJECTS =

utils/liblf.mesh.utils.a: utils/CMakeFiles/lf.mesh.utils.dir/print_info.o
utils/liblf.mesh.utils.a: utils/CMakeFiles/lf.mesh.utils.dir/write_matlab.o
utils/liblf.mesh.utils.a: utils/CMakeFiles/lf.mesh.utils.dir/write_tikz.o
utils/liblf.mesh.utils.a: utils/CMakeFiles/lf.mesh.utils.dir/build.make
utils/liblf.mesh.utils.a: utils/CMakeFiles/lf.mesh.utils.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/hg/u/magina/Documents/lehrfempp/lib/lf/mesh/utils/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Linking CXX static library liblf.mesh.utils.a"
	cd /hg/u/magina/Documents/lehrfempp/lib/lf/mesh/utils/utils && $(CMAKE_COMMAND) -P CMakeFiles/lf.mesh.utils.dir/cmake_clean_target.cmake
	cd /hg/u/magina/Documents/lehrfempp/lib/lf/mesh/utils/utils && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/lf.mesh.utils.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
utils/CMakeFiles/lf.mesh.utils.dir/build: utils/liblf.mesh.utils.a

.PHONY : utils/CMakeFiles/lf.mesh.utils.dir/build

utils/CMakeFiles/lf.mesh.utils.dir/clean:
	cd /hg/u/magina/Documents/lehrfempp/lib/lf/mesh/utils/utils && $(CMAKE_COMMAND) -P CMakeFiles/lf.mesh.utils.dir/cmake_clean.cmake
.PHONY : utils/CMakeFiles/lf.mesh.utils.dir/clean

utils/CMakeFiles/lf.mesh.utils.dir/depend:
	cd /hg/u/magina/Documents/lehrfempp/lib/lf/mesh/utils && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /hg/u/magina/Documents/lehrfempp/lib/lf/mesh /hg/u/magina/Documents/lehrfempp/lib/lf/mesh/utils /hg/u/magina/Documents/lehrfempp/lib/lf/mesh/utils /hg/u/magina/Documents/lehrfempp/lib/lf/mesh/utils/utils /hg/u/magina/Documents/lehrfempp/lib/lf/mesh/utils/utils/CMakeFiles/lf.mesh.utils.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : utils/CMakeFiles/lf.mesh.utils.dir/depend
