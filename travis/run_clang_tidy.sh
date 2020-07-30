# This script is called from .travis.yml
#
# It takes two inputs as environment variables:
#   COMPILER specifies the compiler to use (e.g. g++-7)
#   BUILD_TYPE specifies the CMAKE_BUILD_TYPE.
# 
# It builds and tests Lehrfempp using the provided compiler + cmake configuration.

# Exit immediately from this script upon error
set -e

mkdir -p ${HUNTER_ROOT}

# Initialize cmake and build all dependencies
source $(dirname $0)/build_dependencies.sh

# run clang-tidy
cd Build
cmake -DCMAKE_EXPORT_COMPILE_COMMANDS=On .
../run_clang_tidy.sh