cmake_minimum_required (VERSION 2.8.8)

FIND_PACKAGE(Doxygen)

if(DOXYGEN_FOUND)
  SET(DOXYGEN_INPUT Doxyfile.in)
  SET(DOXYGEN_OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
  SET(DOXYGEN_STRIP_FROM_INC_PATH ${CMAKE_CURRENT_LIST_DIR}/../../lib)

  
  
  OPTION(LF_DOX_EXTRACT_PRIVATE "Equals the value value EXTRACT_PRIVATE of doc/doxygen/doxyfile. Set this to On in order to also generate doxygen documentation of private members." ON)
  OPTION(LF_DOX_INCLUDE_PROJECTS "Set this to true to include all the projects in the ./projects folder in the doxygen documentation." OFF)
  OPTION(LF_DOX_ENABLE_LREF "If this is set to false, @lref{} statements in doxygen documentation are not translated into formatted latex label." ON)
  
  if(LF_DOX_EXTRACT_PRIVATE)
    set(LF_DOX_EXTRACT_PRIVATE_ YES)
  else()
    set(LF_DOX_EXTRACT_PRIVATE_ NO)
  endif()
  
  if(LF_DOX_INCLUDE_PROJECTS)
    set(LF_DOX_INCLUDE_PROJECTS_ "@CMAKE_CURRENT_SOURCE_DIR@/../../projects")
  else()
    set(LF_DOX_INCLUDE_PROJECTS_ "")
  endif()
  
  if(LF_DOX_ENABLE_LREF)
    # Setup filtering of doxygen source files so that we can turn Latex references 
    # from NPDELrefs.aux into section numbers in the doxygen doc.
    
    # C++ script which is executed for every doxygen source file:
    set(CMAKE_CXX_FLAGS_DEBUG ${CMAKE_CXX_FLAGS_RELEASE}) # compile in release mode, even if building debug (for efficiency)
    set(CMAKE_CXX_FLAGS_MINSIZEREL ${CMAKE_CXX_FLAGS_RELEASE})
    add_executable(doxygen_filter filter.cpp)
    target_link_libraries(doxygen_filter PUBLIC Boost::boost) 
    set_target_properties(doxygen_filter PROPERTIES 
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}$<0:hello>"
      RUNTIME_OUTPUT_NAME "doxygen_filter"
      CXX_STANDARD 17
    )
    # Name of the file with LaTeX referencing information
    set(doxygen_input_filter "${CMAKE_CURRENT_BINARY_DIR}/doxygen_filter${CMAKE_EXECUTABLE_SUFFIX} ${CMAKE_CURRENT_LIST_DIR}/NPDEFLrefs.aux")
  else()
    set(doxygen_input_filter "")
  endif()
  
  
  
  # request to configure the file
  configure_file(${DOXYGEN_INPUT} ${DOXYGEN_OUTPUT} @ONLY)
  
  # note the option all which allows to build the docs together with the application
  add_custom_target(doxygen 
    COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUTPUT}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating API Documentation with Doxygen"
    VERBATIM
    )
  if(LF_DOX_ENABLE_LREF)
    add_dependencies(doxygen doxygen_filter)
  endif()

  endif()
  
  add_subdirectory(snippets)
