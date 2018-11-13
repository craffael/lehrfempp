# include <iostream>
# include "lf/base/comm.h"
# include "outside.h"

using namespace lf::base;

void sub() {
  cc::Debug(0, "Sub: A level 0 debug message");
  cc::Debug(10, "Sub: A level 10 debug message");

  std::cout << "Sub: x = " << cv::Get("x") << "\n";
  std::cout << "Sub: y = " << cv::Get("y") << "\n";
}

int main(int argc, char** argv) {
  if (argc <= 1) {
    std::cout << "Demo for lf::base::comm.\n" \
                 " * Try setting the variables v1,v2,v3,v4 by --v1 <val>\n" \
                 " * Set verbosity with the -v option\n" \
                 " * Print debug messages, by --debug_levels or --debug_code, "\
                 "available levels: 0,1,2,5,10. Set levels directly by\n" \
                 "   --debug_levels 1,2,5\n" \
                 "   or via the binary code: 1,2,5=100110=38, hence\n" \
                 "   --debug_code 38\n\n";
  }
  // initialise input reading
  ci::Init(argc, argv); // argc, argv optional. can be specified later

  // add some variables
  ci::Add<int>("v1", "Set variable #1"); // directly
  ci::Add<double>("v2", "Set variable #2", 42.1); // with default value
  ci::Add("help", "Print this help message.");  // add help message

  bool verbose = false;
  ci::Add() // or via the boost::program_options interface
    ("v3", po::value<int>()->default_value(10), "Set variable #3.")
    ("v4", po::value<std::string>(), "Set variable #4.")
    ("verbose,v", po::bool_switch(&verbose), "Set verbosity.");

  // parse command line & file for variables
  ci::ParseCommandLine(); // can also take argc, argv as arguments
  ci::ParseFile("params.par"); // get value from params.par

  ci::Help();  // print help & var names, if -h was specified
  if (verbose)
    std::cout << "Verbosity on.\n";
  else
    std::cout << "Verbose not set.\n";

  // print only if "v1" is set. If it isn't set and Get is called a invalid
  // argument exception is thrown
  if (ci::IsSet("v1"))
    std::cout << "Value of v1 is " << ci::Get<int>("v1") << "\n";

  // Safe, since it has a default value
  std::cout << "Value of v2 is " << ci::Get<double>("v2") << "\n";
  // get "v4", but if it isn't set it return "empty". Used to avoid errors.
  std::cout << "Value of v4 is " << ci::Get<std::string>("v4", "empty") << "\n";
  // std::cout << ci::Get<double>("v5") << "\n"; // invalid_argument exception!

  // verify a certain statement
  cc::Verify(ci::Get<double>("v2") == 42.1, "v2 must be 42.1!");

  // debug messages of different levels
  cc::Debug(0, "A level 0 debug message.");
  cc::Debug(1, "A level 1 debug message.");
  cc::Debug(2, "A level 2 debug message.");

  // global variables
  cv::Add("x", 2);
  cv::Add("y", -100);
  std::cout << "Globally set 'x'=2 and 'y'=-100\n";

  // the global variables are also accesible from othre functions & files
  sub();
  outside();

  return 0;
}
