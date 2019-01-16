# include "comm.h"

namespace lf::base {

namespace comm {

namespace variables {

std::map<std::string, std::pair<bs::hold_any, std::string>> kGlobalVars;

void ListVariables() {
  for (auto const& el : kGlobalVars) {
    std::string key = el.first;
    auto val = el.second; // nested pairs. Clearer than el.second.second.
    std::cout << key << " = " << val.first;
    if (val.second != "") // if description not empty
      std::cout << " : " << val.second;
    std::cout << "\n";
  }
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

  // check for StaticVars defined
  // ctrl_root globally defined in static_vars.cc
  const StaticVar* it = ctrl_root;
  while (it != nullptr) {
    // Avoid duplicate entries
    try {
      // Don't add the option if it exists already (then no error is thrown)
      // false -> only exact match in name is admissible
      po::option_description el = kDesc.find(it->name_, false);
    }
    catch (const std::exception e) {
      Add()
        (it->name_.c_str(), po::value<unsigned int>(&it->ref_), it->comment_.c_str());
    }
    it = it->next_;
  }
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
void ParseCommandLine(const int argc, const char** argv) {
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


} // namespace input

} // namespace comm

} // namepace lf::base
