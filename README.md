[![Build Status](https://travis-ci.org/craffael/lehrfempp.svg?branch=master)](https://travis-ci.org/craffael/lehrfempp)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/craffael/lehrfempp)](https://ci.appveyor.com/project/craffael/lehrfempp)

# LehreFEM++
Simplistic Finite Element Framework for research and eduction optimzed for clarity and flexibility with some trade-off concerning performance.

## Naming conventions
* Filenames should be all lowercase and can include underscores (_) or dashes (-).
* Type names start with a capital letter and have a capital letter for each new word, with no underscores.
* The names of variables (including function parameters) and data members are all lowercase, with underscores between words.
* Data members of classes (not structs), both static and non-static, are named like ordinary nonmember variables, but with a trailing underscore.
* Regular functions have mixed case; accessors and mutators may be named like variables. Ordinarily, functions should start with a capital letter and have a capital letter for each new word.
* Namespace names are all lower-case.
* Macros are all uppercase.
* Variables declared constexpr or const, and whose value is fixed for the duration of the program, are named with a leading "k" followed by mixed case.

[Doxygen Class Documentation](https://craffael.github.io/lehrfempp)
