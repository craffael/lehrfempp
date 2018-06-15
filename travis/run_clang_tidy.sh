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
pwd
source ./install_cmake.sh

# compile
cd ${TRAVIS_BUILD_DIR}
export CXX=clang++-6.0
cmake -H. -BBuild -DCMAKE_BUILD_TYPE=Debug -D -Wdev

# run clang-tidy
$(dirname $0)/run-clang-tidy.py -j2 -p ${TRAVIS_BUILD_DIR} -header-filter=lib/ '^((?!snippets).)*(?<!_tests\.cc)$'