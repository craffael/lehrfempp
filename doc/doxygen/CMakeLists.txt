
FIND_PACKAGE(Doxygen)

if(DOXYGEN_FOUND)
  SET(DOXYGEN_INPUT Doxyfile.in)
  SET(DOXYGEN_OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
  SET(DOXYGEN_STRIP_FROM_INC_PATH ${CMAKE_CURRENT_LIST_DIR}/../../lib)

  
  
  OPTION(LF_DOX_EXTRACT_PRIVATE "Equals the value value EXTRACT_PRIVATE of doc/doxygen/doxyfile. Set this to On in order to also generate doxygen documentation of private members." ON)
  OPTION(LF_DOX_INCLUDE_PROJECTS "Set this to true to include all the projects in the ./projects folder in the doxygen documentation." OFF)
  OPTION(LF_DOX_ENABLE_FILTERS "If this is set to false, filters (see filters/ directory) are not run" ON)

  
  if(LF_DOX_EXTRACT_PRIVATE)
    set(LF_DOX_EXTRACT_PRIVATE_ YES)
  else()
    set(LF_DOX_EXTRACT_PRIVATE_ NO)
  endif()
  
  if(LF_DOX_INCLUDE_PROJECTS)
    set(LF_DOX_INCLUDE_PROJECTS_ "${CMAKE_CURRENT_SOURCE_DIR}/../../projects")
  else()
    set(LF_DOX_INCLUDE_PROJECTS_ "")
  endif()
  
  if(LF_DOX_ENABLE_FILTERS)
    add_subdirectory(filters)
    # Name of the file with LaTeX referencing information
    set(DOXYGEN_FILTER "${CMAKE_CURRENT_BINARY_DIR}/filters/doxygen_filter${CMAKE_EXECUTABLE_SUFFIX} --project-root ${CMAKE_SOURCE_DIR} --aux ${CMAKE_CURRENT_LIST_DIR}/NPDEFLrefs.aux")
    message(STATUS "Doxygen filters enabled." ${DOXYGEN_FILTER})
  else()
    set(DOXYGEN_FILTER "")
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
  set_target_properties(doxygen PROPERTIES FOLDER "doc")

  if(LF_DOX_ENABLE_FILTERS)
    add_dependencies(doxygen doxygen_filter)
  endif()

  
  if(LF_BUILD_DOC)
    add_subdirectory(snippets)
  endif()

endif(DOXYGEN_FOUND)