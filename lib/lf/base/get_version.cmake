# Functionality for lf::base::LehrfemppVersion()
find_package(Git)
# if(Git_FOUND)
  # # retrieve commit sha:
  # execute_process(COMMAND ${GIT_EXECUTABLE} status --porcelain=v2 --branch
                  # WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
                  # RESULT_VARIABLE git_result
                  # OUTPUT_VARIABLE git_out
                 # )
  # if(git_result) 
    # # non-zero return code:
    # message(WARNING "Error while running git status --porcelain=v2 --branch from ${CMAKE_CURRENT_LIST_DIR}. lf::base::LehrfemVersion() will not return any information.")
  # else()
    # # Retrieve branch.oid: sha of the current commit
    # if(git_out MATCHES "# branch.oid ([a-z0-9]+)")
      # set(version_sha ${CMAKE_MATCH_1})
    # else()
      # set(version_sha "NOTFOUND")
    # endif()
    
    # message("version_sha = ${version_sha}")
    
    # # retrieve commit date time:
    # execute_process(COMMAND ${GIT_EXECUTABLE} show --no-patch --no-notes "--pretty=%cI" ${version_sha}
                    # WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
                    # RESULT_VARIABLE git_result
                    # OUTPUT_VARIABLE git_out
                   # )
                 
    # if(git_result)
      # # non-zero return code:
      # message(WARNING "Error while running git show --no-patch --no-notes --pretty='%cI' ${version_sha} from ${CMAKE_CURRENT_LIST_DIR}. lf::base::LehrfemVersion() will not return any information.")
    # else()
      # string (STRIP ${git_out} version_datetime)
    # endif()    
  # endif()
# endif()

if(NOT version_sha OR version_sha STREQUAL "NOTFOUND")
  # In this case we try to retrieve the version sha via hunter:
  message("hello: ${__HUNTER_FINAL_URL_lehrfempp}")
endif()

  message(WARNING git not found, lf::base::GitInfo() will not return any information.)


message("version_sha = ${version_sha}")
message("version_datetime = ${version_datetime}")
message("hunter_lehrfempp_version = ${HUNTER_lehrfempp_VERSION}")