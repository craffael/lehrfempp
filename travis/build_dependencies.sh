# This script is called from .travis.yml
#
# It takes two inputs as environment variables:
#   COMPILER specifies the compiler to use (e.g. g++-7)
#   BUILD_TYPE specifies the CMAKE_BUILD_TYPE.
# 
# It initializes cmake and makes sure that hunter downloads and builds all dependencies

# Exit immediately from this script upon error
set -e

# install new version of cmake:
mkdir -p ${HUNTER_ROOT}

# Install cmake
source $(dirname $0)/install_cmake.sh

# compile
cd ${TRAVIS_BUILD_DIR}
export CXX=${COMPILER}
$CXX --version

cmake -H. -BBuild -DCMAKE_BUILD_TYPE=${BUILD_TYPE} -DCMAKE_CXX_FLAGS=-g0 -DHUNTER_CONFIGURATION_TYPES=${BUILD_TYPE} -Wdev
