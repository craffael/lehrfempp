# include <iostream>
# include <fstream>
# include <string>
# include <boost/program_options.hpp> // compile w/ -lboost_program_options
# include "lf/base/base.h" // for StaticVars class

namespace po = boost::program_options;

void parse_file(const std::string& config_file, po::options_description& desc,
                po::variables_map& vm){
  // @brief Parse config_file for variables specified in desc and store in vm
  // @param config_file: file containing the variables in form var=value
  // @param desc: description of variables, see main()
  // @param vm: variables map, stores the variables with a key

  std::ifstream config_fstream(config_file);
  vm = po::variables_map(); // clear variables
  // extract variables from config file and store in vm.
  // it is possible to also allow variables that haven't been declared yet,
  // by the bool value true. this is useful, since otherwise boost throws an
  // error, if the values haven't been declared!
  po::store(po::parse_config_file(config_fstream, desc, true), vm);
  po::notify(vm);
}

void parse_cmd_line(int argc, char** argv, po::options_description& desc,
                    po::variables_map& vm){
  // @brief Parse command line for variables specified in desc and store in vm
  // @param argc: argc from main(argc, argv)
  // @param argv: argv from main(argc, argv)
  // @param desc: description of variables, see main()
  // @param vm: variables map, stores the variables with a key

  vm = po::variables_map(); // clear variables
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);
}

void print_variables(po::options_description& desc, po::variables_map& vm) {
  // @brief Print the variables in vm or print help message in desc
  // @param desc: description of variables, see main()
  // @param vm: variables map, stores the variables with a key

  // check for variables via the count method:
  if (vm.count("help"))
    std::cout << desc << "\n"; // prints info about the variables
  if (vm.count("ctrl_var"))
    std::cout << "ctrl_var:  " << vm["ctrl_var"].as<int>() << "\n";
  if (vm.count("other_var"))
    std::cout << "other_var: " << vm["other_var"].as<int>() << "\n";
}

class StaticVarsDemoClass {
 public:
  static int ctrl_var_;
  static int other_var_;
};
CLASSCONTROLDECLARE(StaticVarsDemoClass, ctrl_var_, "ctrl_var");
CLASSCONTROLDECLARE(StaticVarsDemoClass, other_var_, "other_var");

void set_static_vars(po::options_description& desc, po::variables_map& vm) {
  // @brief Set the static variables in StaticVarsDemoClass, if they are set
  //        in the variables map vm
  // @param desc: description of variables, see main()
  // @param vm: variables map, stores the variables with a key

  // we could now directly set the static variables
  // (if they are specified in the command line, i.e. count("..") > 0)
  if (vm.count("ctrl_var")) {
    std::cout << "Setting static variable: ctrl_var.\n";
    lf::base::SetCtrlVar("ctrl_var", vm["ctrl_var"].as<int>());
  }
  if (vm.count("other_var")){
    std::cout << "Setting static variable: other_var.\n";
    lf::base::SetCtrlVar("other_var", vm["other_var"].as<int>());
  }
  std::cout << "Value of static variable ctrl_var: "
            << StaticVarsDemoClass::ctrl_var_ << "\n"
            << "Value of static variable other_var: "
            << StaticVarsDemoClass::other_var_ << "\n";
}

int main(int argc, char** argv) {
  std::cout << "Demo program for boost::program_options."
            << "Variables can be set via the command line, run\n"
            << "./a.out --help\n"
            << "for a list of arguments. Alternatively the parameters can be"
            << "set via a file, here we refer to setup.vars.\n";

  // add possible program options
  po::options_description desc("Possible options");
  desc.add_options()
    ("help", "produce help message")
    ("ctrl_var", po::value<int>(), "(int) Set variable #1.")
    ("other_var", po::value<int>(), "(int) Set variable #2.");

  // the variable map stores the variables
  po::variables_map vm;

  // load variables via setup.vars
  parse_file("setup.vars", desc, vm);
  std::cout << "\nValues read in from setup.vars:\n";
  // desc is only given as argument to print the help message, if needed
  print_variables(desc, vm);

  // load variables via command line
  parse_cmd_line(argc, argv, desc, vm);
  std::cout << "\nValues read in from command line:\n";
  print_variables(desc, vm);
  std::cout << "\n";

  // set static variables and print them
  set_static_vars(desc, vm);

  return 0;
}
