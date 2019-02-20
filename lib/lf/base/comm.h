#ifndef __comm_h
#define __comm_h

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/program_options.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/spirit/home/support/detail/hold_any.hpp>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>

#include "static_vars.h"

// namespace structure:
// lf::base
// --------::comm
// --------------::variables (for global variables)
// --------------::input (for reading variables from cmdline & files)

namespace lf::base {

namespace bs = boost::spirit;  // faster, templated version of boost::any
namespace po = boost::program_options;  // keep in lf::base to avoid conflicts

namespace comm {

namespace variables {

// (key, (value, comment)) = (string, (hold_any, string))
extern std::map<std::string, std::pair<bs::hold_any, std::string>> kGlobalVars;

template <typename T>
void Add(const std::string& key, const T& value,
         const std::string& comment = "");

template <typename T>
T Get(const std::string& key);

/**
 * @brief List all global variables as "name" "value" pairs
 */
extern void ListVariables();

/**
 * @brief Checks if key exists in kGlobalVars.
 * @param key The key of the (key, value) pair.
 * @return true, if key exists, false otherwise.
 */
extern bool IsSet(const std::string& key);

/**
 * @brief Removes key and return true if it existed, false otherwise.
 * @param key The key of the (key, value) pair.
 * @return true if key existed, false otherwise.
 */
extern bool Remove(const std::string& key);

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

    throw std::invalid_argument("The key " + key + " couldn't be found.");
  }
}

}  // namespace variables

namespace input {

extern int kArgc;
extern char** kArgv;
extern std::string kConfigFile;
extern po::variables_map kVM;
extern po::options_description kDesc;

/**
 * @brief Initialises parameters to read from command line (argc, argv) and
 *        for reading from a file. Also sets default options debug_code and
 *        debug_levels.
 * @param argc argc from `int main(int argc, char** argv)`.
 * @param argv argv from `int main(int argc, char** argv)`.
 * @param file A file containing options in form name=value.
 */
extern void Init(int argc, char** argv, const std::string& file);

/**
 * @brief Interface to input(argc, argv, file). Sets file=std::string().
 * @param argc: int, argc from `int main(int argc, char** argv)`.
 * @param argv: char**, argv from `int main(int argc, char** argv)`.
 */
extern void Init(int argc, char** argv);

/**
 * @brief Interface to input(argc, argv, file). Sets argc=0, argv=nullptr.
 * @param file A file containing options in form name=value.
 */
extern void Init(const std::string& file);

/**
 * @brief Interface to input(argc, argv, file).
 *        Sets argc=0, argv=nullptr, file=std:string().
 */
extern void Init();

/**
 * @brief Handle to po::options_description.add_options.
 * @return desc.add_options()
 * @note Use as: `ci::add()("v1", po::value<int>(), "Description")(..);
 */
extern po::options_description_easy_init Add();

/**
 * @brief Add (name, comment) pair to options description.
 * @param name The name of variable.
 * @param comment Description of the variable.
 * @note Equivalent to `desc.add_options()(name.c_str(), comment.c_str())`.
 */
extern void Add(const std::string& name, const std::string& comment);

/**
 * @brief Add possible input for variable called `name` with description
 *        `comment`.
 * @param name The name of the variable.
 * @param comment Description of the variable.
 */
template <typename T>
void Add(const std::string& name, const std::string& comment);

/**
 * @brief Add possible input for variable called `name` with description
 *        `comment` with default value `def`.
 * @param name The name of the variable.
 * @param comment Description of the variable.
 * @param def T The default value.
 */
template <typename T>
void Add(const std::string& name, const std::string& comment, const T& def);

/**
 * @brief Possibility of setting any variable `value`. Option is called `name`
          with description `comment`
 * @param name The name of the option.
 * @param value The variable to set.
 * @param comment (optional) Description of the option.
 */
template <typename T>
void AddSetter(const std::string& name, T& value,
               const std::string& comment = "");

/**
 * @brief Get the value of the variable `name`.
 * @param name The name of the variable.
 * @note Throws an invalid_argument exception if `name` doens't exist.
 */
template <typename T>
T Get(const std::string& name);

/**
 * @brief Get the value of the variable `name`, return `alt` if `name` doesn't
 *        exist. 'Safe' version of Get<T>(name).
 * @param name The name of the variable.
 * @param alt T Value that's returned, if `name` doesn't exist.
 */
template <typename T>
T Get(const std::string& name, const T& alt);

/**
 * @brief Print help, if it exists in options_description.
 * @return true, if "help" was set, false otherwise.
 */
extern bool Help();

/**
 * @brief Check if option "name" has been set.
 * @param name The name of variable for which to check.
 * @return true, if found, false otherwise.
 */
extern bool IsSet(const std::string& name);

/**
 * @brief Parse (argc, argv) for values of options.
 * @param argc argc from `int main(int argc, char** argv)`. (optional)
 * @param argv argv from `int main(int argc, char** argv)`. (optional)
 * @note If argc and argv are not provided, the values given to init will
 *       be used. If none have been given, no variables will be found.
 */
extern void ParseCommandLine(const int& argc = 0, const char** argv = nullptr);

/**
 * @brief Parse the config file for variables of form name=value.
 * @param file A file with variables. (optional)
 * @return true if file exists (kConfigFile if file = ""), false if not.
 * @note If file is not provided, the file given to init will be used.
 *       If none has been given, no variables will be found.
 */
extern bool ParseFile(const std::string& file = "");

template <class T>
class Track {
 public:
  /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  Constructor: Places a new item of the global info list in the list
  Usually called via a macro (DECLARE, COUNTER).
  The second version of the constructor also permits to add a comment to the
  information item.
  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  Track(const std::string& name, T& ref,
        const std::string& comment = std::string());
  Track(const std::string& name, T& ref, const T& def,
        const std::string& comment = std::string());
  Track() = delete;
  Track(const Track&) = delete;
  Track(Track&&) = delete;
  Track& operator=(const Track&) = delete;
  Track& operator=(Track&&) = delete;
  ~Track() = default;
};

template <class T>
Track<T>::Track(const std::string& name, T& ref, const std::string& comment) {
  try {
    // Don't add the option if it exists already (then no error is thrown)
    // false -> only exact match in name is admissible
    const po::option_description& el = kDesc.find(name, false);
  } catch (const std::exception& e) {
    Add()(name.c_str(), po::value<unsigned int>(&ref), comment.c_str());
  }
}

template <class T>
Track<T>::Track(const std::string& name, T& ref, const T& def,
                const std::string& comment) {
  try {
    // Don't add the option if it exists already (then no error is thrown)
    // false -> only exact match in name is admissible
    const po::option_description& el = kDesc.find(name, false);
  } catch (const std::exception& e) {
    Add()(name.c_str(), po::value<unsigned int>(&ref)->default_value(def),
          comment.c_str());
  }
}

// Type for managing static variables
using StaticVar = Track<unsigned int>;

template <typename T>
void Add(const std::string& name, const std::string& comment) {
  kDesc.add_options()(name.c_str(), po::value<T>(), comment.c_str());
}

template <typename T>
void Add(const std::string& name, const std::string& comment, const T& def) {
  kDesc.add_options()(name.c_str(), po::value<T>()->default_value(def),
                      comment.c_str());
}

namespace internal {  // functions that are not supposed to be accessed by the
                      // user

/**
 * @brief Returns a lambda function to set the variable given to this function
 * @return See brief.
 */
template <typename T>
std::function<void(T)> SetValue(T& value) {
  std::function<void(T)> lambda = [&value](T new_value) { value = new_value; };
  return lambda;
}

}  // namespace internal

template <typename T>
void AddSetter(const std::string& name, T& value, const std::string& comment) {
  // Avoid duplicate entries
  try {
    // Don't add the option if it exists already (then no error is thrown)
    // false -> only exact match in name is admissible
    const po::option_description& el = kDesc.find(name, false);
  } catch (const std::exception& e) {
    // If not found, then we add the option
    kDesc.add_options()(name.c_str(),
                        po::value<T>()->notifier(internal::SetValue(value)),
                        comment.c_str());
  }
}

template <typename T>
T Get(const std::string& name) {
  if (kVM.count(name) > 0) return kVM[name].as<T>();
  throw std::invalid_argument(
      "In template Get<T>(const std::string&): "
      "Value ``" +
      name + "'' not set. Terminating.");
  return T();
}

template <typename T>
T Get(const std::string& name, const T& alt) {
  try {
    return Get<T>(name);
  } catch (const std::invalid_argument& e) {
    return alt;
  }
}

}  // namespace input

}  // namespace comm

namespace ci = comm::input;
namespace cv = comm::variables;

}  // namespace lf::base

/**
 * @brief Create a new element of type lf::base::Track<unsigned>
 *        with the given variable, name and description.
 *        This will later be used to add a command line option called
 *        "name" for setting the variable "uintvar" with the
 *        description "comment".
 * @param uintvar The variable we can set from command line.
 * @param name What the option will be called (--<name>)
 * @param comment The description of the option.
 */
#define ADDOPTION(uintvar, name, comment) \
  unsigned int uintvar = 0;               \
  static lf::base::ci::Track<unsigned int> name(#name, uintvar, comment)

/**
 * @brief Create a new element of type lf::base::Track<unsigned>
 *        with the given variable, name and description.
 *        This will later be used to add a command line option called
 *        "name" for setting the variable "uintvar" with the
 *        description "comment".
 * @param uintvar The variable we can set from command line.
 * @param default The default value for the variable.
 * @param name What the option will be called (--<name>)
 * @param comment The description of the option.
 */
#define ADDOPTION_DEFAULT(uintvar, default, name, comment)               \
  unsigned int uintvar = 0; /* could also set it =default */             \
  static lf::base::ci::Track<unsigned int> name(#name, uintvar, default, \
                                                comment)

#endif  // __comm_h
