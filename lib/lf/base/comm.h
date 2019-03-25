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

#include "lf_assert.h"
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

/**
 * @brief Add a new global variable triple with (key, value, comment)
 * @param key Key of the variable
 * @param value Value of the variable
 * @param comment Comment to the variable (optional)
 */
template <typename T>
void Add(const std::string& key, const T& value,
         const std::string& comment = "");

/**
 * @brief Retrieve the value of a global variable for a given key
 * @param key Key of the variable
 * @return Value of the variable to the given key
 */
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
  }
  throw std::invalid_argument("The key " + key + " couldn't be found.");
}

}  // namespace variables

namespace input {

extern std::string& getConfigFile();
extern po::variables_map& getVM();
extern po::options_description& getDesc();

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
 * @param argc argc from `int main(int argc, char** argv)`.
 * @param argv argv from `int main(int argc, char** argv)`.
 */
extern void ParseCommandLine(const int& argc, char** argv);

/**
 * @brief Parse the config file for variables of form name=value.
 * @param file A file with variables.
 * @return true if file exists (kConfigFile if file = ""), false if not.
 */
extern bool ParseFile(const std::string& file);

template <class T>
class Track {
 public:
  /**
   * @brief Constructor for (option_name, &value_to_be_changed, comment)
   *        triple
   * @param name The name of the option
   * @param ref If the option is set it will change the value of ref
   * @param comment A description of the option
   */
  Track(const std::string& name, T& ref,
        const std::string& comment = std::string());
  /**
   * @brief Constructor for (option_name, &value_to_be_changed,
   *        default_value, comment) quadrupel
   * @param name The name of the option
   * @param ref If the option is set it will change the value of ref
   * @param ref The default value of the option
   * @param comment A description of the option
   */
  Track(const std::string& name, T& ref, const T& def,
        const std::string& comment = std::string());

  // delete default & copy constructor, not needed
  Track() = delete;
  Track(const Track&) = delete;
  Track(Track&&) = delete;
  Track& operator=(const Track&) = delete;
  Track& operator=(Track&&) = delete;

  /**
   * @brief The default destructor
   */
  ~Track() = default;
};

template <class T>
Track<T>::Track(const std::string& name, T& ref, const std::string& comment) {
  auto temp = getDesc().find_nothrow(name.c_str(), false);

// TODO(craffael): Remove the following line if
// https://developercommunity.visualstudio.com/content/problem/297876/static-inline-variable-gets-destroyed-multiple-tim.html
// is fixed
#ifndef _MSC_VER
  LF_ASSERT_MSG(getDesc().find_nothrow(name.c_str(), false) == nullptr,
                "Name conflict: There is already another boost program option "
                "with the name " +
                    name + " registered.");
#endif

  Add()(name.c_str(), po::value<unsigned int>(&ref), comment.c_str());
}

template <class T>
Track<T>::Track(const std::string& name, T& ref, const T& def,
                const std::string& comment) {
// TODO(craffael): Remove the following line if
// https://developercommunity.visualstudio.com/content/problem/297876/static-inline-variable-gets-destroyed-multiple-tim.html
// is fixed
#ifndef _MSC_VER
  LF_ASSERT_MSG(getDesc().find_nothrow(name.c_str(), false) == nullptr,
                "Name conflict: There is already another boost program option "
                "with the name " +
                    name + " registered.");
#endif
  Add()(name.c_str(), po::value<unsigned int>(&ref)->default_value(def),
        comment.c_str());
}

// Type for managing static variables
using StaticVar = Track<unsigned int>;

template <typename T>
void Add(const std::string& name, const std::string& comment) {
// TODO(craffael): Remove the following line if
// https://developercommunity.visualstudio.com/content/problem/297876/static-inline-variable-gets-destroyed-multiple-tim.html
// is fixed
#ifndef _MSC_VER
  LF_ASSERT_MSG(getDesc().find_nothrow(name.c_str(), false) == nullptr,
                "Name conflict: There is already another boost program option "
                "with the name " +
                    name + " registered.");
#endif
  getDesc().add_options()(name.c_str(), po::value<T>(), comment.c_str());
}

template <typename T>
void Add(const std::string& name, const std::string& comment, const T& def) {
// TODO(craffael): Remove the following line if
// https://developercommunity.visualstudio.com/content/problem/297876/static-inline-variable-gets-destroyed-multiple-tim.html
// is fixed
#ifndef _MSC_VER
  LF_ASSERT_MSG(getDesc().find_nothrow(name.c_str(), false) == nullptr,
                "Name conflict: There is already another boost program option "
                "with the name " +
                    name + " registered.");
#endif
  getDesc().add_options()(name.c_str(), po::value<T>()->default_value(def),
                          comment.c_str());
}

template <typename T>
void AddSetter(const std::string& name, T& value, const std::string& comment) {
// TODO(craffael): Remove the following line if
// https://developercommunity.visualstudio.com/content/problem/297876/static-inline-variable-gets-destroyed-multiple-tim.html
// is fixed
#ifndef _MSC_VER
  LF_ASSERT_MSG(getDesc().find_nothrow(name.c_str(), false) == nullptr,
                "Name conflict: There is already another boost program option "
                "with the name " +
                    name + " registered.");
#endif

  auto SetValue = [&value](T new_value) { value = new_value; };
  getDesc().add_options()(name.c_str(), po::value<T>()->notifier(SetValue),
                          comment.c_str());
}

template <typename T>
T Get(const std::string& name) {
  if (getVM().count(name) > 0) {
    return getVM()[name].as<T>();
  }
  throw std::invalid_argument(
      "In template Get<T>(const std::string&): "
      "Value ``" +
      name + "'' not set. Terminating.");
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
#define ADDOPTION(uintvar, name, comment)                         \
  inline unsigned int uintvar = 0;                                \
  inline lf::base::ci::Track<unsigned int> ____comm_##name##____( \
      #name, uintvar, comment)

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
#define ADDOPTION_DEFAULT(uintvar, default, name, comment)          \
  inline unsigned int uintvar = 0; /* could also set it =default */ \
  inline lf::base::ci::Track<unsigned int> ____comm_##name##____(   \
      #name, uintvar, default, comment)

/**
 * @brief Macro for threshold-conditional output
 * @param ctrlvar integer control variable
 * @param level control level
 * @statement code to be executed
 *
 * The code passed in statement is executed if the value of the
 * control variable is larger than the value passed in level
 *
 * @note The executable code must not involve a comma operator.
 * Commas inside strings are ok.
 */
#define CONTROLLEDSTATEMENT(ctrlvar, level, statement) \
  if ((ctrlvar) >= (level)) {                          \
    statement;                                         \
  }

/**
 * @brief Macro for bit-flag-conditional output
 * @param ctrlvar integer control variable
 * @param flagpat selection bit pattern for flags
 * @statement code to be executed
 *
 * The code passed in statement is executed if the value of the
 * control variable is larger than the value passed in level
 *
 * @note The executable code must not involve a comma operator.
 * Commas inside strings are ok.
 */
#define SWITCHEDSTATEMENT(ctrlvar, flagpat, statement) \
  if (((ctrlvar) & (flagpat)) > 0) {                   \
    statement;                                         \
  }

#endif  // __comm_h
