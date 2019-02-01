/** @file static_vars.h
 *
 * For handling of global control variables
 */

#ifndef __c3c605c9e48646758bf03fab65d5283D
#define __c3c605c9e48646758bf03fab65d5283D

#include <string>
#include <utility>

namespace lf::base {
/*
 * @brief Object in the list of global static variables
 *
 */
template <class T>
class Track {
 public:
  Track *next_;         /**< next variable in list */
  Track *prev_;         /**< previous variable in list */
  std::string name_;    /**< name of variable */
  std::string comment_; /**< optional documentation */
  T &ref_;

  /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  Constructor: Places a new item of the global info list in the list
  Usually called via a macro (DECLARE, COUNTER).
  The second version of the constructor also permits to add a comment to the
  information item.
  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  Track(std::string name, T &ref, Track<T> *&root,
        std::string comment = std::string({}));
  Track() = delete;
  Track(const Track &) = delete;
  Track(Track &&) = delete;
  Track &operator=(const Track &) = delete;
  Track &operator=(Track &&) = delete;
  ~Track();

  static int count;
};

template <class T>
int Track<T>::count = 0;

template <class T>
Track<T>::Track(std::string name, T &ref, Track<T> *&root, std::string comment)
    : prev_(nullptr),
      name_(std::move(name)),
      comment_(std::move(comment)),
      ref_(ref) {
  if (root) {
    root->prev_ = this;
  }
  next_ = root;
  root = this;
  count++;
}

template <class T>
Track<T>::~Track() {
  if (prev_) {
    prev_->next_ = next_;
  }
  if (next_) {
    next_->prev_ = prev_;
  }
  count--;
}

// Type for managing static variables
using StaticVar = Track<unsigned int>;
}  // namespace lf::base

/**
 * @name Macros for handling diagnostics control variables
 */
/**@{*/

#define ADDOPTION(uintvar, name, comment)                             \
  unsigned int uintvar = 0;                                           \
  static lf::base::Track<unsigned int> name(#name, uintvar,            \
                                       lf::base::ctrl_root, comment)

#define CONTROLDECLARE(intvar, varname)                       \
  unsigned int intvar = 0;                                    \
  static lf::base::StaticVar ctrlvar##intvar(varname, intvar, \
                                             lf::base::ctrl_root)
#define CONTROLDECLAREINFO(intvar, varname, info)             \
  unsigned int intvar = 0;                                    \
  static lf::base::StaticVar ctrlvar##intvar(varname, intvar, \
                                             lf::base::ctrl_root, #info)
#define EXTERNDECLAREINFO(intvar, varname, info)              \
  extern unsigned(intvar);                                    \
  static lf::base::StaticVar ctrlvar##intvar(varname, intvar, \
                                             lf::base::ctrl_root, #info)

#define CLASSCONTROLDECLARE(class, intvar, varname)                 \
  unsigned int class ::intvar = 0;                                  \
  static lf::base::StaticVar class##intvar(varname, class ::intvar, \
                                           lf::base::ctrl_root)

#define CONTROLDECLARECOMMENT(class, intvar, varname, comment)      \
  unsigned int class ::intvar = 0;                                  \
  static lf::base::StaticVar class##intvar(varname, class ::intvar, \
                                           lf::base::ctrl_root, comment)
                                      // Why was this #comment and not comment?
/**@}*/

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

namespace lf::base {

/**
 * @brief Pointer to global container for control variables
 *
 * This pointer is defined in lf_assert.cc
 */
extern StaticVar *ctrl_root;

/**
 * @brief Output of globally managed static integers
 *
 */
unsigned int ListCtrlVars(std::ostream &out,
                          const StaticVar *ctrl_var_root = nullptr);

/**
 * @brief Setting a named control variable
 * @param varname string containing the name of the global control variable
 * @param ctrl_var_root pointer to the beginning of the variable list
 */
bool SetCtrlVar(const std::string &varname, int value,
                const StaticVar *ctrl_var_root = nullptr);

/**
 * @brief Initialize variables from file
 *
 * File format is `varname=value` per line, where `value` has to be a string
 * that can be converted into an integer.
 * @param filename name of the file to read from
 */
bool ReadCtrlVarsFile(const std::string &filename,
                      const StaticVar *ctrl_var_root = nullptr);

/** @brief Process command line arguments to set variables
 *
 * Command line arguments of the form
 *    `vaname=value'
 * can be used to set internal static integer variables.
 * The variable name `varname` must not contain the character `=`
 */
int ReadCtrVarsCmdArgs(int argc, const char *argv[],
                       const StaticVar *ctrl_var_root = nullptr);
}  // end namespace lf::base

#endif  // __c3c605c9e48646758bf03fab65d52836
