cmake_minimum_required (VERSION 2.8.8)

FIND_PACKAGE(Doxygen)

if(DOXYGEN_FOUND)
  SET(DOXYGEN_INPUT Doxyfile.in)
  SET(DOXYGEN_OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
  SET(DOXYGEN_STRIP_FROM_INC_PATH ${CMAKE_CURRENT_LIST_DIR}/../../lib)

  # Settings for filtering of doxygen comments
  # Name of filter script (Python3 script)
  SET(filter_script filter.py)
  SET(filter_script_output ${CMAKE_CURRENT_BINARY_DIR}/${filter_script})
  # Name of the file with LaTeX referencing information
  SET(aux_file NPDEFLrefs.aux)
  SET(aux_file_output ${CMAKE_CURRENT_BINARY_DIR}/${aux_file})
  SET(doxygen_filter ${filter_script})
  SET(doxygen_aux ${aux_file})
  configure_file(${filter_script} ${filter_script_output} COPYONLY)
  configure_file(${aux_file} ${aux_file_output} COPYONLY)
  
  OPTION(LF_DOX_EXTRACT_PRIVATE "Equals the value value EXTRACT_PRIVATE of doc/doxygen/doxyfile. Set this to On in order to also generate doxygen documentation of private members." ON)
  OPTION(LF_DOX_INCLUDE_PROJECTS "Set this to true to include all the projects in the ./projects folder in the doxygen documentation." OFF)
  
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
  
  # request to configure the file
  configure_file(${DOXYGEN_INPUT} ${DOXYGEN_OUTPUT} @ONLY)
  
  # note the option all which allows to build the docs together with the application
  add_custom_target(doxygen 
    COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUTPUT}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating API Documentation with Doxygen"
    VERBATIM
    )

  endif()
  
  add_subdirectory(snippets)
