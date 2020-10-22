# Functionality for lf::base::GitInfo()
find_package(Git)
if(Git_FOUND)
  execute_process(COMMAND ${GIT_EXECUTABLE} status --porcelain=v2 --branch
                  WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
                  RESULT_VARIABLE git_result
                  OUTPUT_VARIABLE git_out
                 )
  if(git_result) 
    # non-zero return code:
    message(WARNING "Error while running git status --porcelain=v2 --branch from ${CMAKE_CURRENT_LIST_DIR}. lf::base::GitInfo() will not return any information.")
  else()
  message("#\\sbranch.oid\\s([0-9a-z]{40})")
    # Retrieve branch.oid: sha of the current commit
    if(git_out MATCHES "# branch.oid ([a-z0-9]+)")
      set(git_branch_oid ${CMAKE_MATCH_1})
    else()
      set(git_branch_oid "NOTFOUND")
    endif()
    
    # Retrieve branch.head: The name of the current branch
    if(git_out MATCHES "# branch.head ([0-9a-z\\/\\-\\.]+)")
      set(git_branch_head ${CMAKE_MATCH_1})
    else()
      set(git_branch_head "NOTFOUND")
    endif()
    
    # Retrieve branch.upstream: The name of the upstream branch (if set)
    if(git_out MATCHES "# branch.upstream ([0-9a-z\\/\\-\\.]+)")
      set(git_branch_upstream ${CMAKE_MATCH_1})
    else()
      set(git_branch_upstream "NOTFOUND")
    endif()
    
    message("git_branch_oid = ${git_branch_oid}")
    message("git_branch_head = ${git_branch_head}")
    message("git_branch_upstream = ${git_branch_upstream}")
    
  endif()
  
else()
  message(WARNING git not found, lf::base::GitInfo() will not return any information.)
endif()

message("hunter_lehrfempp_version = ${HUNTER_lehrfempp_VERSION}")