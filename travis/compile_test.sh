# This script is called from .travis.yml
#
# It takes two inputs as environment variables:
#   COMPILER specifies the compiler to use (e.g. g++-7)
#   BUILD_TYPE specifies the CMAKE_BUILD_TYPE.
# It builds and tests Lehrfempp using the provided compiler + cmake configuration.

# Exit immediately from this script upon error
set -e

# Initialize cmake and build all dependencies
source $(dirname $0)/build_dependencies.sh

cd Build
make -j${NUM_PROC:-2}

# test
ctest --output-on-failure