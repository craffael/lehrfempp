# include "comm.h"

namespace lf::base {

namespace comm {

namespace variables {

std::map<std::string, int> kGlobalVars;

/**
 * @brief Add (key, value) pair to kGlobalVars map.
 * @param key The key of the (key, value) pair.
 * @param value The value of the (key, value) pair.
 */
void Add(const std::string& key, const int value) {
  kGlobalVars[key] = value;
}

/**
 * @brief Get value to corresponding to key.
 * @param key The key of the (key, value) pair.
 * @return The value of the (key, value) pair.
 * @note: Throws error if key is not found. Use with cv::IsSet to avoid
 *        errors: if cv::IsSet(key) { val = cv::Get(key); }
 */
int Get(const std::string& key) {
  auto it = kGlobalVars.find(key);
  if (it != kGlobalVars.end()) {
    return it->second;
  }
  throw "The key " + key + " couldn't be found.";
  return 0; // against compiler warnings
}

/**
 * @brief Checks if key exists in kGlobalVars.
 * @param key The key of the (key, value) pair.
 * @return true, if key exists, false otherwise.
 */
bool IsSet(const std::string& key) {
  return kGlobalVars.find(key) != kGlobalVars.end();
}

/**
 * @brief Removes key and return true if it existed, false otherwise.
 * @param key The key of the (key, value) pair.
 * @return true if key existed, false otherwise.
 */
bool Remove(const std::string& key) {
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

/**
 * @brief Counts necessary bits to represent the integer `i`.
 * @param i The number for which the min. number of bits are counted.
 * return Min. number of necessary bits to represent `i`.
 */
int CountBits(const int i) {
  int n_bits = 1;
  for (int j = 1; j < i; n_bits++) {
    j += 1 << n_bits;
  }
  return n_bits;
}

/**
 * @brief Given an input string levels, set all bits of the return value to
 *        one at the positions specified in levels.
 * @param levels The positions of bits that will be set to 1.
 * @note Use ":" as range operator. Open ends and beginnings allowed.
 *       Count starts at 0 and all ranges are incl. the first & last value.
 *       Examples: * "0,1,3" --> 1011 --> 11
 *                 * "2:4" --> 11100
 *                 * ":" --> all ones, std::numeric_limits<int>::max()
 *                 * ":3" --> 1111
 */
int ExtractDebugCode(const std::string& levels) {
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

/**
 * @brief Print message if the corresponding debug level was set.
 * @param level The debug level of message.
 * @param message The debug message.
 */
void Debug(const int level, const std::string& message) {
  const int internal_level = 1 << level; // 2^level
  if ((kDebugCode & internal_level) > 0) {
    std::cout << message << "\n";
  }
}

/**
 * @brief Set the debug code as integer.
 * @param debug_code The debug code.
 * @note if debug levels are 0,2,5 then debug_code = binary(100101) = 37.
 */
void SetDebugCode(const int debug_code) {
  kDebugCode = debug_code;
}

/**
 * @brief Return current debug code.
 * @return The debug code in integer representation.
 */
int GetDebugCode() {
  return kDebugCode;
}

/**
 * @brief Verify that condition is true. If it's not, throw a runtime error
 *        with `message`.
 * @param condition If condition is false, the error is thrown.
 * @param message The message of the runtime error.
 */
void Verify(const bool condition, const std::string& message) {
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

/**
 * @brief Initialises parameters to read from command line (argc, argv) and
 *        for reading from a file. Also sets default options debug_code and
 *        debug_levels.
 * @param argc argc from `int main(int argc, char** argv)`.
 * @param argv argv from `int main(int argc, char** argv)`.
 * @param file A file containing options in form name=value.
 */
void Init(int argc, char** argv, std::string file) {
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

/**
 * @brief Interface to input(argc, argv, file). Sets file=std::string().
 * @param argc: int, argc from `int main(int argc, char** argv)`.
 * @param argv: char**, argv from `int main(int argc, char** argv)`.
 */
void Init(int argc, char** argv) {
  Init(argc, argv, std::string());
}

/**
 * @brief Interface to input(argc, argv, file). Sets argc=0, argv=nullptr.
 * @param file A file containing options in form name=value.
 */
void Init(std::string file) {
  Init(0, nullptr, file);
}

/**
 * @brief Interface to input(argc, argv, file).
 *        Sets argc=0, argv=nullptr, file=std:string().
 */
void Init() {
  Init(0, nullptr, std::string());
}

/**
 * @brief Handle to po::options_description.add_options.
 * @return desc.add_options()
 * @note Use as: `ci::add()("v1", po::value<int>(), "Description")(..);
 */
extern po::options_description_easy_init Add() {
  return kDesc.add_options();
}

/**
 * @brief Add (name, comment) pair to options description.
 * @param name The name of variable.
 * @param comment Description of the variable.
 * @note Equivalent to `desc.add_options()(name.c_str(), comment.c_str())`.
 */
void Add(const std::string& name, const std::string& comment) {
  kDesc.add_options()(name.c_str(), comment.c_str());
}

/**
 * @brief Print help, if it exists in options_description.
 * @return true, if "help" was set, false otherwise.
 */
bool Help() {
  if (kVM.count("help") > 0) {
    std::cout << kDesc << "\n";
    return true;
  }
  return false;
}

/**
 * @brief Check if option "name" has been set.
 * @param name The name of variable for which to check.
 * @return true, if found, false otherwise.
 */
bool IsSet(const std::string& name) {
  return kVM.count(name) > 0;
}

/**
 * @brief Parse (argc, argv) for values of options.
 * @param argc argc from `int main(int argc, char** argv)`. (optional)
 * @param argv argv from `int main(int argc, char** argv)`. (optional)
 * @note If argc and argv are not provided, the values given to init will
 *       be used. If none have been given, no variables will be found.
 */
void ParseCommandLine(int argc, char** argv) {
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

/**
 * @brief Parse the config file for variables of form name=value.
 * @param file A file with variables. (optional)
 * @note If file is not provided, the file given to init will be used.
 *       If none has been given, no variables will be found.
 */
void ParseFile(const std::string& file) {
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

/**
 * @brief Make a int-variable global in the sense of cv::Add
 * @param name The name of the variable, as specified in ci::Add
 * @return true, if the variable existed and is int, false otherwise.
 */
bool MakeGlobal(const std::string& name) {
  if (IsSet(name)) {
    int val;
    try {
      val = Get<int>(name);
    }
    catch (const boost::bad_any_cast& e) {
      // Get calls vm[name].as<int>, which throws a bad_any_cast
      // error, if the types mismatch (hence no int)
      return false;
    }
    cv::Add(name, val); // int version of Add
  }
  else {
    std::cout << "Variable " << name << " is not set, cannot make it global.\n";
    return false;
  }
}

} // namespace input

} // namespace comm

} // namepace lf::base
