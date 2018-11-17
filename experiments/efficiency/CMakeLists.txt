# cmake configuration file for compilation of LehrFEM++ examples: RefEl

set(runtime_test runtime_test.cc)

add_executable(experiments.efficiency.runtime_test ${runtime_test})

target_link_libraries(experiments.efficiency.runtime_test
  PUBLIC Eigen3::Eigen Boost::boost Boost::timer Boost::chrono Boost::system base)

target_compile_features(experiments.efficiency.runtime_test PUBLIC cxx_std_17)

  
    

  

