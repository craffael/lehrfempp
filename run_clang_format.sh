#!/bin/bash

# This script is a utility script for linux/osx users to run clang-format over all relevant files
# The continuous integration build environment (github actions) runs exactly the same command to check
# if all files are formatted correctly.
#
# HOW TO USE:
# - run this script to reformat all relevant files in this repository.



# find clang-format
if [ -x "$(command -v clang-format)" ]; then
  ct=$(command -v clang-format)
fi
if [ -x "$(command -v clang-format-5.0)" ]; then
  ct=$(command -v clang-format-5.0)
fi
if [ -x "$(command -v clang-format-6.0)" ]; then
  ct=$(command -v clang-format-6.0)
fi
if [ -x "$(command -v clang-format-7)" ]; then
  ct=$(command -v clang-format-7)
fi
if [ -x "$(command -v clang-format-8)" ]; then
  ct=$(command -v clang-format-8)
fi

if [ -z "$ct" ]; then
  echo "clang-format, clang-format-8, clang-format-7, clang-format-6.0 or clang-format-5.0 not found in path"
fi
echo "using $ct"
$(dirname $0)/ci/run-clang-format.py -r --clang-format-executable $ct --color always $(dirname $0)/lib $(dirname $0)/projects
