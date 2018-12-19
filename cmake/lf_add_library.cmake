# \brief  Setup up a new LehrFEM++ module/library and add it to the build-system.
#         Should be called from CMakeLists.txt of a LehrFEMpp module.
# 
# USAGE:
# `lf_add_library(<libname> [HEADER_ONLY] [source1 [source2 ...]] )`
#  
# This adds a new library target (using add_library() ) with the name
# `libname` that is composed of the sources `source1`, `source2`...
# In addition it makes sure that this library target is "exported" 
# so that it can be used by thid party libraries.
#
# \param libname The name of the library target. If the target lies in folder
#                `lib/lf/a/b/` then libname must be `a.b` (everything lowercase)
# \param sourcex A list of header/source files that belong to this target.
#                (similar to how add_library works)
# \param HEADER_ONLY Specify this if `<libname>` is a header only library,
#                    i.e. it must not be compiled.
#
#
# Calling lf_add_library does a few things at once:
# - It creates a new library target (through add_library() ) named `libname` with
#   the specified sources
# - Makes sure that c++17 is used
# - Setup the install location for target `<libname>` by calling install()
# 
#
# \note If your library depends on other cmake targets, you can call
#       `target_link_libraries(<libname> ...)`, target_compile_options(<libname>...)`.
#
# For a simple example of how to use `lf_add_library()` see the file
# `/lib/lf/base/CMakeLists.txt`
#
function(lf_add_library _args)

  # Extract arguments:
  set(_libname ${ARGV0})
  set(_header_only OFF)
  set(_sources ${ARGV}) # The other libraries on which this library depends.
  list(REMOVE_AT _sources 0)
  if(${ARGC} GREATER 1)
    if(${ARGV1} STREQUAL "HEADER_ONLY")
      set(_header_only ON)
      list(REMOVE_AT _sources 0)
    endif()
  endif()
  
  # CHECK INPUT ARGUMENTS:
  #########################################################################
  
  # Check that the library name is lowercase:
  string(TOLOWER ${_libname} _libname_LC)
  if(NOT(${_libname} STREQUAL ${_libname_LC}))
    message(FATAL_ERROR "Error adding Lehrfem++ library '${_libname}': Libname is not lowercase (${_libname})")
  endif()
  
  # Check that the library name starts with lf.a.b.name if this file lies in lib/lf/a/b/
  get_relative_path(_relative_dir ${CMAKE_CURRENT_LIST_DIR})
  if(NOT _relative_dir)
    message(FATAL_ERROR "Error adding LehrFEM++ library '${_libname}': lf_add_library() must be called from a subdirectory of the 'lib/lf' directory.")
  endif()
  string(SUBSTRING ${_relative_dir} 4 -1 _relative_dir)
  
  string(REPLACE "/" "." _libname_test ${_relative_dir})
  if(NOT ${_libname} STREQUAL ${_libname_test})
    message("relative_dir = ${_relative_dir}")
    message(FATAL_ERROR "Error adding LehrFEM++ library '${_libname}': Libname doesn't match with relative directory (${_relative_dir}), libname should be '${_libname_test}'")
  endif()
  
  # Determine base folder of the lf repository:
  string(REPLACE "/lib/${_relative_dir}" "" _path2lf ${CMAKE_CURRENT_LIST_DIR})
  
  # Check that if this is a header-only library no source files are present:
  list(LENGTH ${_libname}_src _srcLen)
  if(${_header_only} AND (${_srcLen} GREATER 0))
    message(FATAL_ERROR "Error adding LehrFEM++ library '${_libname}': This library is HEADER_ONLY but contains sources.")
  endif()
  
  
  
  # Add the CMake library:
  #########################################################################
  if(${_header_only}) 
    #add_library(${_libname} INTERFACE)
    ## Add the following target so the sources are displayed in visual studio:
    #add_custom_target(${_libname}. SOURCES ${${_libname}_h})
    #set(_scope INTERFACE)
    
    add_library(${_libname} ${sources} ${_path2lf}/cmake/dummy.cpp)
    set(_scope PUBLIC)
  else()
    add_library(${_libname} ${sources})
    set(_scope PUBLIC)
  endif()
  target_include_directories(${_libname} PUBLIC
    "$<BUILD_INTERFACE:${LOCAL_INCLUDE_DIRECTORY}>"
    "$<BUILD_INTERFACE:${LOCAL_GENERATED_INCLUDE_DIRECTORY}>"
    $<INSTALL_INTERFACE:include>
  )
  target_compile_features(${_libname} PUBLIC cxx_std_17)
  set_target_properties(${_libname} PROPERTIES PREFIX "${CMAKE_STATIC_LIBRARY_PREFIX}")
  # set_target_properties(${_libname} PROPERTIES FOLDER "Libraries")
  

  install(TARGETS ${_libname} EXPORT LFTargets
    LIBRARY DESTINATION "lib"
    ARCHIVE DESTINATION "lib")

  # Filter out the header file from ${_sources} so we can move them to the include folder  
  set(_headers)
  foreach(s ${_sources})
    if(s MATCHES "(.h|.hpp|.inc)$")
      set(_headers ${_headers} ${s})
    endif()
  endforeach()
  install(FILES ${_headers} DESTINATION "include/${_relative_dir}")
  
  # Register this library:
  # set(LF_ALL_TARGETS ${_libname} ${LF_ALL_TARGETS} CACHE INTERNAL "" FORCE)
endfunction()