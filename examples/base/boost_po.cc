# include <iostream>
# include <fstream>
# include <string>
# include <boost/program_options.hpp>
// compile w/ -lboost_program_options

namespace po = boost::program_options;

void parse_file(const std::string& config_file, po::options_description& desc,
                po::variables_map& vm){
  std::ifstream config_fstream(config_file);
  vm = po::variables_map(); // clear variables
  // extract variables from config file and store in vm
  po::store(po::parse_config_file(config_fstream, desc), vm);
  po::notify(vm);
}

void parse_cmd_line(int argc, char** argv, po::options_description& desc,
                    po::variables_map& vm){
  vm = po::variables_map(); // clear variables
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);
}

void print_variables(po::options_description& desc, po::variables_map& vm) {
  // check for variables via the count method:
  if (vm.count("help"))
    std::cout << desc << "\n";
  if (vm.count("id"))
    std::cout << "id:    " << vm["id"].as<int>() << "\n";
  if (vm.count("name"))
    std::cout << "name:  " << vm["name"].as<std::string>() << "\n";
  if (vm.count("score"))
    std::cout << "score: " << vm["score"].as<double>() << "\n";
}

int main(int argc, char** argv) {
  // add possible program options
  po::options_description desc("Possible options");
  desc.add_options()
    ("help", "produce help message")
    ("id", po::value<int>(), "(int) Set variable #1.")
    ("name", po::value<std::string>(), "(std::string) Set variable #2.")
    ("score", po::value<double>(), "(double) Set variable #3.");

  // the variable map stores the variables
  po::variables_map vm;

  // load variables via setup.conf
  parse_file("setup.conf", desc, vm);
  std::cout << "\nValues read in from setup.conf:\n";
  // desc is only given as argument to print the help message, if needed
  print_variables(desc, vm);

  // load variables via command line
  parse_cmd_line(argc, argv, desc, vm);
  std::cout << "\nValues read in from command line:\n";
  print_variables(desc, vm);

  return 0;
}
