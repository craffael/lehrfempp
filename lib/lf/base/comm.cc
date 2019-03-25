#include "comm.h"

namespace lf::base {

comm::input::StaticVar* ctrl_root = nullptr;

namespace comm {

namespace variables {

std::map<std::string, std::pair<bs::hold_any, std::string>> kGlobalVars;

void ListVariables() {
  for (auto const& el : kGlobalVars) {
    std::string key = el.first;
    auto val = el.second;  // nested pairs. Clearer than el.second.second.
    std::cout << key << " = " << val.first;
    if (!val.second.empty()) {  // if description not empty
      std::cout << " : " << val.second;
    }
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

std::string& getConfigFile() {
  static std::string value;
  return value;
}
po::variables_map& getVM() {
  static po::variables_map value;
  return value;
};
po::options_description& getDesc() {
  static po::options_description value("Allowed options");
  return value;
}

extern po::options_description_easy_init Add() {
  return getDesc().add_options();
}

void Add(const std::string& name, const std::string& comment) {
  getDesc().add_options()(name.c_str(), comment.c_str());
}

bool Help() {
  if (getVM().count("help") > 0) {
    std::cout << getDesc() << "\n";
    return true;
  }
  return false;
}

bool IsSet(const std::string& name) { return getVM().count(name) > 0; }

void ParseCommandLine(const int& argc, char** argv) {
  po::store(po::parse_command_line(argc, argv, getDesc()), getVM());
  po::notify(getVM());
}

bool ParseFile(const std::string& file) {
  // if file is set, use it. otherwise use this->kConfigFile
  std::ifstream config_fs(!file.empty() ? file : getConfigFile());
  if (!config_fs.good()) {
    return false;
  }

  po::store(po::parse_config_file(config_fs, getDesc()), getVM());
  po::notify(getVM());
  return true;
}

}  // namespace input

}  // namespace comm

}  // namespace lf::base
