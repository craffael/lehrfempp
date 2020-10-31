# This script is called from .travis.yml
#
# It takes two inputs as environment variables:
#   COMPILER specifies the compiler to use (e.g. g++-7)
#   BUILD_TYPE specifies the CMAKE_BUILD_TYPE.
# 
# Also it takes exactly one command line parameter which is a regex expression
# that specifies the files to include (see .travis.yml for example)
# It builds and tests Lehrfempp using the provided compiler + cmake configuration.

# Exit immediately from this script upon error
set -e

mkdir -p ${HUNTER_ROOT}

# Initialize cmake and build all dependencies
source $(dirname $0)/build_dependencies.sh

# run clang-tidy
cd Build
cmake -DCMAKE_EXPORT_COMPILE_COMMANDS=On .
# Run clang-tidy
#../run_clang_tidy.sh --files $1
clang-tidy-8 -p . ../projects/ipdg_stokes/examples/lid_driven_cavity/vortex.cc