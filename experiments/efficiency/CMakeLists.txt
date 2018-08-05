# cmake configuration file for compilation of LehrFEM++ examples: RefEl

set(runtime_test runtime_test.cc)

add_executable(lf.experiments.efficiency.runtime_test ${runtime_test})

target_link_libraries(lf.experiments.efficiency.runtime_test
  PUBLIC Eigen3::Eigen Boost::boost lf.base)

target_compile_features(lf.experiments.efficiency.runtime_test PUBLIC cxx_std_17)

  
    

  

