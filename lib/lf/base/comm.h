# ifndef __comm_h
# define __comm_h

# include <iostream>
# include <fstream>
# include <sstream>
# include <map>
# include <stdexcept>
# include <boost/spirit/home/support/detail/hold_any.hpp>
# include <boost/program_options.hpp>
# include <boost/program_options/options_description.hpp>
# include <boost/algorithm/string/split.hpp>
# include <boost/algorithm/string/classification.hpp>

# include "lf/base/static_vars.h"


// namespace structure:
// lf::base
// --------::comm
// --------------::variables (for global variables)
// --------------::input (for reading variables from cmdline & files)


namespace lf::base {

namespace bs = boost::spirit; // faster, templated version of boost::any
namespace po = boost::program_options; // keep in lf::base to avoid conflicts

namespace comm {

namespace variables {

// (key, (value, comment)) = (string, (hold_any, string))
extern std::map<std::string, std::pair<bs::hold_any, std::string>> kGlobalVars;

template <typename T>
void Add(const std::string&, const T&, const std::string& comment = "");

template <typename T>
T Get(const std::string&);

extern void ListVariables();

extern bool IsSet(const std::string&);

extern bool Remove(const std::string&);

/**
 * @brief Add a the value `value` to the global variables with key `key` Handle
 *        description `comment`.
 * @param key The key of the variable.
 * @param value The value.
 * @param comment (optional) Description of the variable.
 */
template <typename T>
void Add(const std::string& key, const T& value, const std::string& comment) {
  kGlobalVars[key] = std::make_pair(bs::hold_any(value), comment);
}

/**
 * @brief Get the value of the variable `key`.
 * @param key The key of the variable.
 * @note Throws an invalid_argument exception if `key` doens't exist.
 */
template <typename T>
T Get(const std::string& key) {
  if (kGlobalVars.count(key)) {
    return bs::any_cast<T>(kGlobalVars[key].first);
  }
  else {
    throw std::invalid_argument("The key " + key + " couldn't be found.");
  }
}

} // namespace variables

namespace input {

extern int kArgc;
extern char** kArgv;
extern std::string kConfigFile;
extern po::variables_map kVM;
extern po::options_description kDesc;

extern void Init(int, char**, std::string);
extern void Init(int, char**);
extern void Init(std::string);
extern void Init();

extern po::options_description_easy_init Add();
extern void Add(const std::string&, const std::string&);
template <typename T>
void Add(const std::string&, const std::string&);
template <typename T>
void Add(const std::string&, const std::string&, const T&);
template <class A>
void AddCtrl(const std::string&, A&, const std::string& comment = "");
template <typename T>
void AddSetter(const std::string&, T&, const std::string& comment = "");

template <typename T>
T Get(const std::string&);
template <typename T>
T Get(const std::string&, const T&);

extern bool Help();

extern bool IsSet(const std::string&);

extern void ParseCommandLine(const int argc = 0, const char** argv = nullptr);

extern bool ParseFile(const std::string& file = "");

/**
 * @brief Add possible input for variable called `name` with description
 *        `comment`.
 * @param name The name of the variable.
 * @param comment Description of the variable.
 */
template <typename T>
void Add(const std::string& name, const std::string& comment) {
  kDesc.add_options()(name.c_str(), po::value<T>(), comment.c_str());
}

/**
 * @brief Add possible input for variable called `name` with description
 *        `comment` with default value `def`.
 * @param name The name of the variable.
 * @param comment Description of the variable.
 * @param def T The default value.
 */
template <typename T>
void Add(const std::string& name, const std::string& comment, const T& def) {
  kDesc.add_options()
    (name.c_str(), po::value<T>()->default_value(def), comment.c_str());
}

template <class A>
std::function<void(unsigned int)> SetCtrl(A& a) {
  std::function<void(unsigned int)> lambda = [&a](unsigned int c){ a.output_ctrl_ = c; };
  return lambda;
}

/**
 * @brief Possibility of setting the public member `ctrl` of the
 *        `class_instance`. Option is called `name` with description `comment`
 * @param name The name of the option.
 * @param class_instance Instance of the class wher to change member value.
 * @param comment (optional) Description of the option.
 */
template <class A>
void AddCtrl(const std::string& name, A& class_instance,
             const std::string& comment) {
  // Why not solve it like:
  //   void AddCtrl(const string& name, const T& val, const string& comment);
  // Call it like:
  //   AddCtrl("ctrl", MyClass.PublicVar, "Some comment");
  //   AddCtrl("ctrl", MyClass::StaticVar, "Some comment");
  // ?
  kDesc.add_options()
    (name.c_str(),
     po::value<unsigned int>()->notifier(SetCtrl(class_instance)),
     comment.c_str());
}


template <typename T>
std::function<void(T)> SetValue(T& value) {
  std::function<void(T)> lambda = [&value](T new_value){ value = new_value; };
  return lambda;
}

/**
 * @brief Possibility of setting any variable `value`. Option is called `name`
          with description `comment`
 * @param name The name of the option.
 * @param value The variable to set.
 * @param comment (optional) Description of the option.
 */
template <typename T>
void AddSetter(const std::string& name, T& value, const std::string& comment) {
  // Avoid duplicate entries
  try {
    // Don't add the option if it exists already (then no error is thrown)
    // false -> only exact match in name is admissible
    po::option_description el = kDesc.find(name, false);
  }
  catch (const std::exception e) {
    // If not found, then we add the option
    kDesc.add_options()
      (name.c_str(), po::value<T>()->notifier(SetValue(value)), comment.c_str());
  }
}

/**
 * @brief Get the value of the variable `name`.
 * @param name The name of the variable.
 * @note Throws an invalid_argument exception if `name` doens't exist.
 */
template <typename T>
T Get(const std::string& name) {
  if (kVM.count(name) > 0) {
    return kVM[name].as<T>();
  }
  else {
    throw std::invalid_argument("In template Get<T>(const std::string&): " \
                                "Value ``" + name + "'' not set. Terminating.");
    return T();
  }
}

/**
 * @brief Get the value of the variable `name`, return `alt` if `name` doesn't
 *        exist. 'Safe' version of Get<T>(name).
 * @param name The name of the variable.
 * @param alt T Value that's returned, if `name` doesn't exist.
 */
template <typename T>
T Get(const std::string& name, const T& alt) {
  try {
    return Get<T>(name);
  }
  catch (const std::invalid_argument& e) {
    return alt;
  }
}

} // namespace input

} // namespace comm

namespace ci = comm::input;
namespace cv = comm::variables;

} // namespace lf::base

# endif // __comm_h
