/** @file static_vars.h
 *
 * For handling of global control variables
 */

#ifndef __c3c605c9e48646758bf03fab65d5283D
#define __c3c605c9e48646758bf03fab65d5283D

#include <string>
#include <utility>

/**
 * @name Macros for handling diagnostics control variables
 */
/**@{*/

/*
#define CONTROLDECLARE(intvar, varname)                       \
  unsigned int intvar = 0;                                    \
  static lf::base::StaticVar ctrlvar##intvar(varname, intvar, \
                                             lf::base::ctrl_root)
#define CONTROLDECLAREINFO(intvar, varname, info)             \
  unsigned int intvar = 0;                                    \
  static lf::base::StaticVar ctrlvar##intvar(varname, intvar, \
                                             lf::base::ctrl_root, #info)
#define EXTERNDECLAREINFO(intvar, varname, info)              \
  extern unsigned intvar;                                     \
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
*/
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

#endif  // __c3c605c9e48646758bf03fab65d52836
