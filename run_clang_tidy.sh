#!/bin/bash

# This script is a utility script for linux/osx users to run clang-tidy over all relevant files
# The continuous integration build environment (github actions) runs exactly the same command.
#
# HOW TO USE:
# - change into your CMake build directory and call this script from there.
# - If you want to use a custom regex for filtering the files on which clang-tidy is run, add the `--files` option.

set -e

# Check that compile_commands.json exists.
if [[ ! -e ./compile_commands.json ]]; then
  echo "compile_commands.json not found in current directory"
  echo "run 'cmake -DCMAKE_EXPORT_COMPILE_COMMANDS=On . ' to create it."
  exit 1
fi

# find clang-tidy
if [ -x "$(command -v clang-tidy)" ]; then
  ct=clang-tidy
fi
if [ -x "$(command -v clang-tidy-17)" ]; then
  ct=clang-tidy-17
fi

if [ -z "$ct" ]; then
  echo "clang-tidy or clang-tidy-17 not found in path"
  exit 1
fi
echo $ct
version=$($ct --version)
if [[ ! $version =~ "LLVM version 17." ]]; then
  echo "ERROR: Found clang-tidy but it doesn't have version 17.x"
  exit 1
fi

# Parse commandline arguments to see if --files option has been passed:
FILES='^((?!snippets|/test/|/test_utils/).)*$'
POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -f|--files)
    FILES="$2"
    shift # past argument
    shift # past value
    ;;
    *)    # unknown option
    POSITIONAL+=("$1") # save it in an array for later
    shift # past argument
    ;;
esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters

$(dirname $0)/ci/run-clang-tidy.py -p . -clang-tidy-binary $ct -header-filter=lib/ $FILES $1
