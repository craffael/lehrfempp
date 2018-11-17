
# \brief Takes any path in which a component of the form "lib/lf" appears and extracts
#        just this relative component. E.g. "/home/raffael/KTI/lehrfem/lehrfem/lib/lf/base/fem.h" is transformed
#        to "lib/lf/base/base.h"
#
# \param out [out]   Specify the name of variable into which the relative path should be written.
#                    If `path` does not contain "lib/lf/" this variable will be set to FALSE.
# \param path [in]   The full path from which the relative component should be extracted.
#
# Example usage:
# \code
# get_relative_path(rel "/home/raffael/KTI/lehrfem/lib/lf/base/base.h")
# message(${rel}) # will return lib/lf/base/base.h
# \endcode
#
#
function(get_relative_path out path)
  
  # Match regex to extract anything that starts with lib/lf/
  string(REGEX MATCH "\/lib\/lf[A-Za-z0-9\/_]*" _relative_dir ${path})
  if(NOT _relative_dir)
    set(${out} FALSE PARENT_SCOPE)
	return()
  endif()
  
  #message("input = ${path}")
  # remove all lib/lf/ occurences until only one is left:
  while(TRUE)
    string(SUBSTRING ${_relative_dir} 7, -1 _relative_dir)
	  if(NOT _relative_dir)
	    break() #Also break if _relative_dir is an empty string
	  endif()
	  string(REGEX MATCH "\/lib\/lf(\/[A-Za-z0-9\/_]*)?" _temp ${_relative_dir})
	  if(NOT _temp)
	    break()
	  endif()
	  set(_relative_dir ${_temp})
  endwhile(TRUE)
  if("${_relative_dir}" STREQUAL "")
    set(${out} "lib/lf" PARENT_SCOPE)
	#message("out = hydi, relativedir = ${_relative_dir}")
  else()
    set(${out} "lib/lf${_relative_dir}" PARENT_SCOPE)
	#message("out = hydi${_relative_dir}")
  endif()
  
endfunction()