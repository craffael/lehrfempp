set(sources
    intergridfe.h
    mesh_function_transfer.h
)

lf_add_library(lf.intergridfe ${sources})
target_link_libraries(lf.intergridfe PUBLIC Eigen3::Eigen lf.base)

if (LF_ENABLE_TESTING)
    add_subdirectory(test)
endif()
