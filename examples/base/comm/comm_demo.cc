# include <iostream>
# include "lf/base/base.h"
# include "outside.h"

using namespace lf::base;

// a class with a static member which we will set from 
// the command line 
class Dummy {
public:
  // We will set ctrl from the command line later.
  static unsigned int ctrl;
};
unsigned int Dummy::ctrl; // define the variable

// Accessing global variables from another function
// than in which they were added (which is main)
void sub() {
  std::cout << "Sub: x = " << cv::Get<int>("x") << "\n";
  std::cout << "Sub: y = " << cv::Get<std::string>("y") << "\n";
}

// must be char**, not const char** or (const) char* []
int main(int argc, char** argv) {
  // namespace overview:
  // ci (comm::input) for reading in variables from cmdline or file
  // cv (comm::variables) for global variables

  if (argc <= 1) {
    std::cout << "Demo for lf::base::comm. Try the help option "\
                 "(-h, --help) for more instructions.\n";
  }
  // initialise input reading
  ci::Init(argc, argv); // argc, argv optional, can also be specified later

  // - Adding variables
  //   To read in from the command line (or file) we add variables 
  //   to comm::input (or simply `ci`)
  // -- Add help message, this should always be done!
  ci::Add("help,h", "Print this help message.");  
  // -- Add some variables via:
  //    ci::Add<T>("name", "description", [optional] default value)
  ci::Add<int>("v1", "Set variable #1", 0); // with default value
  ci::Add<double>("v2", "Set variable #2"); // no default value
  // -- Or via the boost::program_options interface:
  bool verbose = false;
  ci::Add()
    ("v3", po::value<int>()->default_value(10), "Set variable #3.")
    ("v4", po::value<std::string>(), "Set variable #4.")
    ("verbose,v", po::bool_switch(&verbose), "Set verbosity.");
 
  // - Add the possibility to set existing variables from 
  //   the commend line (or file):
  //     ci::AddSetter<T>("name", variable, "comment")
  //   So the variable verbose could also be set this way:
  //     ci::AddSetter<bool>("verbose,v", verbose, "Set verbosity.");
  //   Set the member variable ctrl of Dummy
  ci::AddSetter<unsigned int>("ctrl", Dummy::ctrl, "Set static member ctrl class Dummy.");

  // - Parse command line & file for variables
  ci::ParseCommandLine(); // can also take argc, argv as arguments
  ci::ParseFile("params.par"); // get value from params.par

  // - Check for the set variables
  // -- Is the help option set?
  if (ci::Help())  // print help & var names, if -h was specified
    return 0;      // and stop here! No need to go through the program :)

  // -- Check if verbose has been set to true
  if (verbose)
    std::cout << "Verbosity on.\n";
  else
    std::cout << "Verbose not set.\n";
 
  // -- Print the value of "v2" if it has been set.
  //    If it isn't set and Get is called a invalid argument exception is thrown
  if (ci::IsSet("v2"))
    std::cout << "Value of v2 is " << ci::Get<double>("v2") << "\n";
  else
    std::cout << "Variable v2 is not set.\n";

  // -- Print the value of v1. This is safe, since it has a default value
  std::cout << "Value of v1 is " << ci::Get<int>("v1") << "\n";
  // -- Get "v4", but if it isn't set it return "empty". Used to avoid errors.
  std::cout << "Value of v4 is " << ci::Get<std::string>("v4", "empty") << "\n";
  // std::cout << ci::Get<double>("v5") << "\n"; // invalid_argument exception!

  // -- Print value of Dummy::ctrl
  std::cout << "Value of Dummy::ctrl is " << Dummy::ctrl<< "\n";

  // - Global variables
  // -- Add an integer and a string
  cv::Add("x", 2);
  cv::Add("y", std::string("Applepie."));
  std::cout << "Globally set x=2 and y='Applepie.'\n";

  // -- Check for the variables
  //    The global variables are also accesible from other functions & files.
  sub();
  outside();

  return 0;
}
