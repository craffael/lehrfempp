# Functionality for lf::base::LehrFemInfo:
include(cmake/get_version.cmake)
include(CheckCXXSourceCompiles)
configure_file(lehrfem_info.cc ${CMAKE_CURRENT_BINARY_DIR}/lehrfem_info.cc)

# Define sources that make up this Library:
set(sources
  base.h
  eigen_tools.h
  invalid_type_exception.h
  lehrfem_info.h
  ${CMAKE_CURRENT_BINARY_DIR}/lehrfem_info.cc
  lf_assert.cc
  lf_assert.h
  lf_exception.h
  predicate_true.h
  ref_el.cc
  ref_el.h
  scalar_concept.h
  spdlog_utils.h
  spdlog_utils.cc
  timer.h
  timer.cc
  types.h
)

lf_add_library(lf.base ${sources})
target_link_libraries(lf.base PUBLIC Eigen3::Eigen Boost::boost spdlog::spdlog)
set_target_properties(lf.base PROPERTIES FOLDER "lib")

# required libraries for boost::stacktrace:
if(UNIX)
  # check if we can compile lf_assert.cc with backtrace support:
  check_cxx_source_compiles(${CMAKE_CURRENT_LIST_DIR}/cmake/has_backtrace.cc assert_compiles)
  
  if(NOT assert_compiles)
    # We have to find the backtrace.h file ourselves:
    if(NOT BACKTRACE_H)
      # we have to supply the path to the backtrace include separately
      # see https://www.boost.org/doc/libs/1_75_0/doc/html/stacktrace/configuration_and_build.html
      # 1) if gcc is installed, lets try to get the include paths from there:
      execute_process(
                      COMMAND /usr/bin/gcc -E -Wp,-v -xc /dev/null
                      RESULT_VARIABLE GCC_ERROR
                      ERROR_VARIABLE GCC_OUTPUT
                      OUTPUT_QUIET
                      )
      set(HINTS "")
      if(NOT GCC_ERROR)
        # extract include paths using regex:
        if(GCC_OUTPUT MATCHES "#include <\\.\\.\\.> search starts here:(.*)End of search list")
          separate_arguments(HINTS UNIX_COMMAND ${CMAKE_MATCH_1})
        endif()
      endif()
      
      # Try to find backtrace.h
      find_file(BACKTRACE_H "backtrace.h" 
        HINTS
        ${HINTS}
      )
    endif()
    if(NOT BACKTRACE_H)
      message(STATUS "Your compiler can't find backtrace.h => backtraces will not include source lines.")
      target_link_libraries(lf.base PUBLIC dl)
    else()
      message(STATUS "Found backtrace.h - ${BACKTRACE_H}")
      target_compile_definitions(lf.base PRIVATE "BOOST_STACKTRACE_BACKTRACE_INCLUDE_FILE=\"${BACKTRACE_H}\"")
      target_compile_definitions(lf.base PRIVATE BOOST_STACKTRACE_USE_BACKTRACE)
      target_link_libraries(lf.base PUBLIC dl backtrace)
    endif()
    if(APPLE)
      # The following is additionally needed on OSX
      target_compile_definitions(lf.base PRIVATE _GNU_SOURCE)
    endif()
  else()
    # backtrace.h can be included as is:
    target_compile_definitions(lf.base PRIVATE BOOST_STACKTRACE_USE_BACKTRACE)
    target_link_libraries(lf.base PUBLIC dl backtrace)
  endif()
  
endif()

if(LF_REDIRECT_ASSERTS)
  # redirect boost asserts to LF_ASSERT/LF_VERIFY:
  target_compile_definitions(lf.base PUBLIC BOOST_ENABLE_ASSERT_DEBUG_HANDLER)
  target_compile_definitions(lf.base PUBLIC "LF_REDIRECT_ASSERTS")
endif()

# enable SPDLOG_DEBUG and SPDLOG_TRACE statements in all modes except release::
target_compile_definitions(lf.base PUBLIC $<IF:$<CONFIG:Release>,SPDLOG_ACTIVE_LEVEL=SPDLOG_LEVEL_INFO,SPDLOG_ACTIVE_LEVEL=SPDLOG_LEVEL_TRACE>)

if(MSVC)
  if(${MSVC_VERSION} GREATER_EQUAL 1915) 
    # You must acknowledge that you understand MSVC resolved a byte alignment issue in this compiler
    # We get this due to using Eigen objects and allocating those objects with make_shared
    target_compile_definitions( lf.base PUBLIC _ENABLE_EXTENDED_ALIGNED_STORAGE )
  endif()
endif()

if(LF_ENABLE_TESTING)
  add_subdirectory(test)
endif()


