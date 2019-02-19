#include "comm.h"

namespace lf::base {

namespace comm {

namespace variables {

std::map<std::string, std::pair<bs::hold_any, std::string>> kGlobalVars;

void ListVariables() {
  for (auto const& el : kGlobalVars) {
    std::string key = el.first;
    auto val = el.second;  // nested pairs. Clearer than el.second.second.
    std::cout << key << " = " << val.first;
    if (val.second != "")  // if description not empty
      std::cout << " : " << val.second;
    std::cout << "\n";
  }
}

bool IsSet(const std::string& key) {
  return kGlobalVars.find(key) != kGlobalVars.end();
}

bool Remove(const std::string& key) {
  auto it = kGlobalVars.find(key);
  if (it != kGlobalVars.end()) {
    kGlobalVars.erase(it);
    return true;
  }
  return false;
}

}  // namespace variables

namespace input {

int kArgc;
char** kArgv;
std::string kConfigFile;
po::variables_map kVM;
po::options_description kDesc("Allowed options");

void Init(int argc, char** argv, const std::string& file) {
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
      const po::option_description& el = kDesc.find(it->name_, false);
    } catch (const std::exception& e) {
      Add()(it->name_.c_str(), po::value<unsigned int>(&it->ref_),
            it->comment_.c_str());
    }
    it = it->next_;
  }
}

void Init(int argc, char** argv) { Init(argc, argv, std::string()); }

void Init(const std::string& file) { Init(0, nullptr, file); }

void Init() { Init(0, nullptr, std::string()); }

extern po::options_description_easy_init Add() { return kDesc.add_options(); }

void Add(const std::string& name, const std::string& comment) {
  kDesc.add_options()(name.c_str(), comment.c_str());
}

bool Help() {
  if (kVM.count("help") > 0) {
    std::cout << kDesc << "\n";
    return true;
  }
  return false;
}

bool IsSet(const std::string& name) { return kVM.count(name) > 0; }

void ParseCommandLine(const int& argc, const char** argv) {
  // maybe argc/argv haven't been set yet and are given as arguments to this
  if (argc != 0 && argv != nullptr) {
    po::store(po::parse_command_line(argc, argv, kDesc), kVM);
  } else {
    // if they have already been set via constructor, use them
    po::store(po::parse_command_line(kArgc, kArgv, kDesc), kVM);
  }
  po::notify(kVM);
}

bool ParseFile(const std::string& file) {
  // if file is set, use it. otherwise use this->kConfigFile
  std::ifstream config_fs(!file.empty() ? file : kConfigFile);
  if (!config_fs.good()) 
    return false;

  po::store(po::parse_config_file(config_fs, kDesc), kVM);
  po::notify(kVM);
  return true;
}

}  // namespace input

}  // namespace comm

}  // namespace lf::base
