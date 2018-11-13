# include "comm.h"

namespace lf::base {

namespace comm {

namespace variables {

std::map<std::string, int> kGlobalVars;

void Add(const std::string& key, const int value) {
  // @brief Add (key, value) pair to kGlobalVars map
  // @param key: std::string, key of pair
  // @param value: int, value of pair

  kGlobalVars[key] = value;
}

int Get(const std::string& key) {
  // @brief Get value to corresponding to key
  // @param key: std::string, key of pair
  // @return: int, value of pair
  // @note: Throws error if key is not found. Use with cv::IsSet to avoid
  //        errors: if cv::IsSet(key) { val = cv::Get(key); }

  auto it = kGlobalVars.find(key);
  if (it != kGlobalVars.end()) {
    return it->second;
  }
  throw "The key " + key + " couldn't be found.";
  return 0; // against compiler warnings
}

bool IsSet(const std::string& key) {
  // @brief Checks if key exists in kGlobalVars
  // @param key: std::string, key of pair
  // @return: true, if key exists, false otherwise

  return kGlobalVars.find(key) != kGlobalVars.end();
}

bool Remove(const std::string& key) {
  // @brief Removes key and return true if it existed, false otherwise.
  // @param key: Key in std::map
  // @return: true if key existed, false otherwise

  auto it = kGlobalVars.find(key);
  if (it != kGlobalVars.end()) {
    kGlobalVars.erase(it);
    return true;
  }
  return false;
}

} // namespace variables

namespace check {

int kDebugCode = 0;

int CountBits(const int i) {
  // @brief Counts necessary bits to represent the integer i
  // @param i: int, number for which the min. number of bits are counted
  // return: min. number of necessary bits to represent i

  int n_bits = 1;
  for (int j = 1; j < i; n_bits++) {
    j += 1 << n_bits;
  }
  return n_bits;
}

int ExtractDebugCode(const std::string& levels) {
  // @brief Given an input string levels, set all bits of the return value to
  //        one at the positions specified in levels
  // @param levels: positions of bits that will be set to 1.
  // @note: Use ":" as range operator. Open ends and beginnings allowed.
  //        Count starts at 0 and all ranges are incl. the first & last value.
  //        Examples: * "0,1,3" --> 1011 --> 11
  //                  * "2:4" --> 11100
  //                  * ":" --> all ones, std::numeric_limits<int>::max()
  //                  * ":3" --> 1111

  int code = 0;
  const int max_code = std::numeric_limits<int>::max(); // all ones
  char delim; // for dumping the delimiters
  int level;
  int last_level = 0; // needed for ranges (:)
  bool is_range = false; // needed for ranges (:)

  std::stringstream ss(levels);

  // if the first element is the range operator, we need to know!
  if(ss.peek() == ':') {
    is_range = true;
    code = 1 << 0; // add level 0, later we add all from last_level + 1 on
    ss.get(); // discard it
  }

  bool done = false;
  while (!done) {
    // check for int
    if(ss >> level) {
      // if we had a ":" character, which represents a range
      if (is_range) {
        // add levels: last_level to level
        // last_level + 1, since last_level was already added
        for (int i = last_level + 1; i <= level; ++i) {
          code += 1 << i;
        }
        is_range = false;
      }
      else {
        code += 1 << level;
      }
      last_level = level;
    }
    else {
      done = true;
    }
    // check for delimiter
    if(ss >> delim) {
      if (delim == ':') {
        is_range = true;
      } // else nothing, just continue
    }
    else {
      done = true;
    }
  }

  // input: "x:" --> set everything from x to max_bits
  if (is_range && last_level > 0) {
    for (int i = last_level + 1; i < CountBits(max_code); ++i) {
      code += 1 << i;
    }
  }

  // input: ":" --> set all levels
  if (is_range && last_level == 0) {
    code = max_code; // max_code has all bits = 1 (or should at least..)
  }

  return code;
}

void Debug(const int level, const std::string& message) {
  // @brief Print message if the corresponding debug level was set
  // @param level: int, debug level of message
  // @param message: std::string, debug message

  const int internal_level = 1 << level; // 2^level
  if ((kDebugCode & internal_level) > 0) {
    std::cout << message << "\n";
  }
}

void SetDebugCode(const int debug_code) {
  // @brief Set the debug code as integer
  // @param debug_code: int, debug code
  // @note: if debug levels are 0,2,5 then debug_code = binary(100101) = 37

  kDebugCode = debug_code;
}

int GetDebugCode() {
  // @brief Return current debug code
  // return: debug code in integer representation

  return kDebugCode;
}

void Verify(const bool condition, const std::string& message) {
  // @brief Verify that condition is true. If it's not, throw a runtime error
  //        with `message`
  // @param condition: bool, it false, the error is thrown
  // @param message: std::string, message of error

  if (!condition) {
    throw std::runtime_error(message);
  }
}

} // namespace check

namespace input {

int kArgc;
char** kArgv;
std::string kConfigFile;
po::variables_map kVM;
po::options_description kDesc("Allowed options");

void Init(int argc, char** argv, std::string file) {
  // @brief Initialises parameters to read from command line (argc, argv) and
  //        for reading from a file. Also sets default options debug_code and
  //        debug_levels.
  // @param argc: int, argc from `int main(int argc, char** argv)`
  // @param argv: char**, argv from `int main(int argc, char** argv)`
  // @param file: std::string, file containing options in form name=value

  kArgc = argc;
  kArgv = argv;
  kConfigFile = file;
  // add debug option per default
  Add<int>("debug_code,d",
    "Debug code. Output all debug statements of level l, if binary"            \
    "representation of <debug_code> is 1 at bit l."                            \
    "Example: `--debug_code 6` will print all debug messages of levels 2 and 3"\
    ", since 5 in binary is 110.\n"                                            \
    "Don't use in combination with --debug_levels!");

  Add<std::string>("debug_levels,l",
    "Debug levels. Set individual levels l using comma-separated values." \
    "Example: `--debug_levels \"0,2,5\"` will print all debug messages of" \
    "levels 0, 2 and 5.\n Don't use in combination with --debug_code!");
}

void Init(int argc, char** argv) {
  // @brief Interface to input(argc, argv, file). Sets file=std::string().
  // @param argc: int, argc from `int main(int argc, char** argv)`
  // @param argv: char**, argv from `int main(int argc, char** argv)`

  Init(argc, argv, std::string());
}

void Init(std::string file) {
  // @brief Interface to input(argc, argv, file). Sets argc=0, argv=nullptr.
  // @param file: std::string, file containing options in form name=value

  Init(0, nullptr, file);
}

void Init() {
  // @brief Interface to input(argc, argv, file).
  //        Sets argc=0, argv=nullptr, file=std:string().

  Init(0, nullptr, std::string());
}

extern po::options_description_easy_init Add() {
  // @brief Handle to po::options_description.add_options
  // @return: desc.add_options()
  // @note: Use as: `ci::add()("v1", po::value<int>(), "Description")(..);

  return kDesc.add_options();
}

void Add(const std::string& name, const std::string& comment) {
  // @brief Add (name, comment) pair to options description.
  // @param name: std::string, name of variable
  // @param comment: std::string, description of variable
  // @note: Equivalent to `desc.add_options()(name.c_str(), comment.c_str())`

  kDesc.add_options()(name.c_str(), comment.c_str());
}

void Help() {
  // @brief Print help, if it exists in options_description

  if (kVM.count("help") > 0) {
    std::cout << kDesc << "\n";
  }
}

bool IsSet(const std::string& name) {
  // @brief Check if option "name" has been set.
  // @param name: std::string, name of variable for which to check
  // @return: true, if found, false otherwise.

  return kVM.count(name) > 0;
}

void ParseCommandLine(int argc, char** argv) {
  // @brief Parse (argc, argv) for values of options
  // @param argc: int (optional), argc from `int main(int argc, char** argv)`
  // @param argv: char** (optional), argv from `int main(int argc, char** argv)`
  // @note: If argc and argv are not provided, the values given to init will
  //        be used. If none have been given, no variables will be found.

  // maybe argc/argv haven't been set yet and are given as arguments to this
  if (argc != 0 && argv != nullptr) {
    po::store(po::parse_command_line(argc, argv, kDesc), kVM);
  }
  else {
  // if they have already been set via constructor, use them
    po::store(po::parse_command_line(kArgc, kArgv, kDesc), kVM);
  }
  po::notify(kVM);
  // extract & set debug level, if set
  if (IsSet("debug_levels")) {
    int code = cc::ExtractDebugCode(Get<std::string>("debug_levels"));
    cc::SetDebugCode(code);
  }
  else if (IsSet("debug_code")) {
    cc::SetDebugCode(Get<int>("debug_code"));
  }
}

void ParseFile(const std::string& file) {
  // @brief Parse the config file for variables of form name=value
  // @param file: std::string (optional), file with variables
  // @note: If file is not provided, the file given to init will be used.
  //        If none has been given, no variables will be found.

  // if file is set, use it. otherwise use this->kConfigFile
  std::ifstream config_fs(file != "" ? file : kConfigFile);
  po::store(po::parse_config_file(config_fs, kDesc), kVM);
  po::notify(kVM);
  // extract & set debug level, if set
  if (IsSet("debug_levels")) {
    int code = cc::ExtractDebugCode(Get<std::string>("debug_levels"));
    cc::SetDebugCode(code);
  }
  else if (IsSet("debug_code")) {
    cc::SetDebugCode(Get<int>("debug_code"));
  }
}

} // namespace input

} // namespace comm

} // namepace lf::base
