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
if [ -x "$(command -v clang-format-10)" ]; then
  ct=$(command -v clang-format-10)
fi

if [ -z "$ct" ]; then
  echo "clang-format or clang-format-10 not found in path"
  exit 1
fi
version=$($ct --version)
if [[ ! $version =~ "version 10." ]]; then
  echo "ERROR: Found clang-format but it doesn't have version 10.x. Please install clang-format-10"
  exit 1
fi
echo "using $ct"
$(dirname $0)/ci/run-clang-format.py -r --clang-format-executable $ct --color always $(dirname $0)/lib $(dirname $0)/projects
