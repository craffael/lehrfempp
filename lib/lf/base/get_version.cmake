# Functionality for lf::base::LehrfemppVersion()
#
# The code in this file tries to determine three pieces of information about the current code:
# - version_sha -> the git sha value that identifies the (last) commit
# - version_datetime -> The datetime of the (last) commit
# - version_tag -> A git tag that is possibly associated with this commit
#                  (hard to determine, is most of the times not set)
#
# There are two ways in which we try to retrieve this information:
# 1) We check if this is a git repository and retrieve the information using ``git status`.
#    This can fail if the code is directly downloaded from github and not checked out as git repository.
#    This happens most notably when lehrfempp is used by a 3rd party project through hunter.
# 2) If 1) fails, we check if the variable __HUNTER_FINAL_URL_lehrfempp is set (which is the case if LehrFEM++
#    is included in another project via hunter). We assume that __HUNTER_FINAL_URL_lehrfempp contains either
#    - a commit sha, or
#    - a git tag
#    in both cases, we ask github via their api for further information to e.g. retrieve the version_datetime. 



# @brief ask github for the datetime of a commit (identified by its sha)
function(get_commit_datetime sha out)
  set(filename ${CMAKE_CURRENT_BINARY_DIR}/download)
  
  file(DOWNLOAD "https://api.github.com/repos/craffael/lehrfempp/git/commits/${sha}" ${filename})
  file(READ ${filename} response)
  
  if(response MATCHES "\"date\" ?: ?\"([0-9TZ:\\-]+)\"")
    set(${out} ${CMAKE_MATCH_1} PARENT_SCOPE)
  endif()

endfunction()


# @brief retrieve commit sha1 from a git tag
function(get_commit_sha_from_tag tag out)
  set(filename ${CMAKE_CURRENT_BINARY_DIR}/download)
  
  file(DOWNLOAD "https://api.github.com/repos/craffael/lehrfempp/git/refs/tags/${tag}" ${filename})
  file(READ ${filename} response)
  
  if(response MATCHES "\"tag\".*\"url\" ?: ?\"([^\"]+)\"")
    file(DOWNLOAD ${CMAKE_MATCH_1} ${filename})
    file(READ ${filename} response)

    if(response MATCHES "\"url\" ?: ?\"[^\"]*/git/commits/([0-9a-z]+)\"")
      set(${out} ${CMAKE_MATCH_1} PARENT_SCOPE)
    endif()
  endif()

endfunction()

# 1) Check if we can retrieve information using git directly
find_package(Git)
if(Git_FOUND)
  # retrieve commit sha:
  execute_process(COMMAND ${GIT_EXECUTABLE} status --porcelain=v2 --branch
                  WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
                  RESULT_VARIABLE git_result
                  OUTPUT_VARIABLE git_out
                 )
  if(git_result) 
    # non-zero return code:
    message(WARNING "Error while running git status --porcelain=v2 --branch from ${CMAKE_CURRENT_LIST_DIR}. lf::base::LehrfemVersion() will not return any information.")
  else()
    # Retrieve branch.oid: sha of the current commit
    if(git_out MATCHES "# branch.oid ([a-z0-9]+)")
      set(version_sha ${CMAKE_MATCH_1})
    endif()
    
    # retrieve commit date time:
    execute_process(COMMAND ${GIT_EXECUTABLE} show --no-patch --no-notes "--pretty=%cI" ${version_sha}
                    WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
                    RESULT_VARIABLE git_result
                    OUTPUT_VARIABLE git_out
                   )
                 
    if(git_result)
      # non-zero return code:
      message(WARNING "Error while running git show --no-patch --no-notes --pretty='%cI' ${version_sha} from ${CMAKE_CURRENT_LIST_DIR}. lf::base::LehrfemVersion() will not return any information.")
    else()
      string (STRIP ${git_out} version_datetime)
    endif()    
  endif()
endif()

# 2) If the code is not a git repository, check if __HUNTER_FINAL_URL_lehrfempp
#    is set and retrieve missing information via github api
if(NOT version_sha)
  # Debug:
  #set(__HUNTER_FINAL_URL_lehrfempp "https://github.com/craffael/lehrfempp/archive/8bc5a1304261d5f5f3bb74a5808b6e998cb681e6.tar.gz")
  #set(__HUNTER_FINAL_URL_lehrfempp "https://github.com/craffael/lehrfempp/archive/release-0.7.21.tar.gz")
  
  # First we check if the value of __HUNTER_FINAL_URL_lehrfempp is still the same
  if(__lf_hunter_final_url STREQUAL __HUNTER_FINAL_URL_lehrfempp)
    # Cache is still valid:
    set(version_sha ${__lf_version_sha})
    set(version_datetime ${__lf_version_datetime})
    set(version_tag ${__lf_version_tag})
    
  else()
    
    
    # assume that the url contains a commit sha:
    if(__HUNTER_FINAL_URL_lehrfempp MATCHES "lehrfempp\\/archive\\/(.+)\\.tar\\.gz")
      get_commit_datetime(${CMAKE_MATCH_1} version_datetime)
      if(version_datetime)
        set(version_sha ${CMAKE_MATCH_1})
      else()
        # url doesn't contain a commit sha, but hopefully it is a tag?
        get_commit_sha_from_tag(${CMAKE_MATCH_1} version_sha)
        if(version_sha)
          get_commit_datetime(${version_sha} version_datetime)
          set(version_tag ${CMAKE_MATCH_1})
        endif()
      endif()
    endif()
    
    # Set cache variables:
    set(__lf_hunter_final_url ${__HUNTER_FINAL_URL_lehrfempp} CACHE INTERNAL "")
    set(__lf_version_sha ${version_sha} CACHE INTERNAL "")
    set(__lf_version_datetime ${version_datetime} CACHE INTERNAL "")
    set(__lf_version_tag ${version_datetime} CACHE INTERNAL "")
  endif()
  
endif()

if(NOT version_sha)
  message(WARNING "Could not determine version of this code. lf::base::LehrFemVersion() will not return any information.")
endif()


# Debug:
#message("version_sha = ${version_sha}")
# message("version_datetime = ${version_datetime}")
#message("version_tag = ${version_tag}")