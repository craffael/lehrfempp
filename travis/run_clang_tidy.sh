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

# Install cmake
echo pwd
source $(dirname $0)/install_cmake.sh

# compile
cd ${TRAVIS_BUILD_DIR}
export CXX=clang++-6.0
cmake -H. -BBuild -DCMAKE_BUILD_TYPE=Debug -DCMAKE_EXPORT_COMPILE_COMMANDS=On -Wdev

# run clang-tidy
cd Build
../run_clang_tidy.sh