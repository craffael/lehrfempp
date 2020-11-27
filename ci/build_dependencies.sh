# This script is called from github actions
#
# It takes two inputs as environment variables:
#   COMPILER specifies the compiler to use (e.g. g++-7)
#   BUILD_TYPE specifies the CMAKE_BUILD_TYPE.
# 
# It initializes cmake and makes sure that hunter downloads and builds all dependencies

# Exit immediately from this script upon error
set -e

HUNTER_ROOT=${GITHUB_WORKSPACE}/hunter

# install new version of cmake:
mkdir -p ${HUNTER_ROOT}

# Install cmake
source $(dirname $0)/install_cmake.sh

# compile
cd ${GITHUB_WORKSPACE}
export CXX=${COMPILER}
$CXX --version

cmake -H. -BBuild -DCMAKE_BUILD_TYPE=${BUILD_TYPE} -DCMAKE_CXX_FLAGS=-g0 -DHUNTER_CONFIGURATION_TYPES=${BUILD_TYPE} -Wdev
