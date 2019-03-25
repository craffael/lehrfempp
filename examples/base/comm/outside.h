/**
 * @file
 * @brief illustrates usage of base::comm framework.
 * @author Raffael Casagrande
 * @date   2019-03-24 06:09:09
 * @copyright MIT License
 */

#ifndef __55220edea03b4be7a2d856ac0d9d28e0
#define __55220edea03b4be7a2d856ac0d9d28e0

#include <lf/base/base.h>

void outside();

ADDOPTION(toPrint, toPrint,
          "the value that the function outside() prints to the command line.");

#endif  // __55220edea03b4be7a2d856ac0d9d28e0
