# include <iostream>
# include "lf/base/comm.h"
# include "outside.h"

using namespace lf::base;

class Dummy {
public:
  unsigned int output_ctrl_; // We will set ctrl from the command line later.
  static unsigned int ctrl;
};

void sub() {
  std::cout << "Sub: x = " << cv::Get<int>("x") << "\n";
  std::cout << "Sub: y = " << cv::Get<std::string>("y") << "\n";
}

unsigned int Dummy::ctrl;
int main(int argc, char** argv) {
  // namespace overview:
  // ci (comm::input) for reading in variables from cmdline or file
  // cv (comm::variables) for global variables


  if (argc <= 1) {
    std::cout << "Demo for lf::base::comm.\n" \
                 " * Try setting the variables v1,v2,v3,v4 by --v1 <val>\n" \
                 " * Set verbosity with the -v option\n";
  }
  // initialise input reading
  ci::Init(argc, argv); // argc, argv optional, can also be specified later

  // Add some variables, either via these functions..
  ci::Add<int>("v1", "Set variable #1", 0); // with default value
  ci::Add<double>("v2", "Set variable #2"); // no default value
  ci::Add("help,h", "Print this help message.");  // add help message
  // ..or via the boost::program_options interface
  bool verbose = false;
  ci::Add()
    ("v3", po::value<int>()->default_value(10), "Set variable #3.")
    ("v4", po::value<std::string>(), "Set variable #4.")
    ("verbose,v", po::bool_switch(&verbose), "Set verbosity.");

  Dummy dummy;
  // Set the member variable ctrl of instance dummy
  ci::AddCtrl("ctrl", dummy, "Set member ctrl of instance dummy.");

  // parse command line & file for variables
  ci::ParseCommandLine(); // can also take argc, argv as arguments
  ci::ParseFile("params.par"); // get value from params.par

  if (ci::Help())  // print help & var names, if -h was specified
    return 0;      // and stop here! No need to go through the program :)

  if (verbose)
    std::cout << "Verbosity on.\n";
  else
    std::cout << "Verbose not set.\n";

  // print only if "v2" is set. If it isn't set and Get is called a invalid
  // argument exception is thrown
  if (ci::IsSet("v2"))
    std::cout << "Value of v2 is " << ci::Get<double>("v2") << "\n";
  else
    std::cout << "Variable v2 is not set.\n";

  // safe, since it has a default value
  std::cout << "Value of v1 is " << ci::Get<int>("v1") << "\n";
  // get "v4", but if it isn't set it return "empty". Used to avoid errors.
  std::cout << "Value of v4 is " << ci::Get<std::string>("v4", "empty") << "\n";
  // std::cout << ci::Get<double>("v5") << "\n"; // invalid_argument exception!

  // print value of dummy.ctrl
  std::cout << "Value of dummy.output_ctrl_ is " << dummy.output_ctrl_<< "\n";

  // global variables
  cv::Add("x", 2);
  cv::Add("y", std::string("Applepie."));
  std::cout << "Globally set x=2 and y='Applepie.'\n";

  // the global variables are also accesible from other functions & files
  sub();
  outside();

  return 0;
}
