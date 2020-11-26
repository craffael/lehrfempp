# This script installs an up-2-date version of cmake on the travis build server.

# Exit immediately from this script upon error
set -e

# save current directory
dir=$(pwd)

mkdir -p ${DEPS_DIR} && cd ${DEPS_DIR}

if [ ! -d "doxygen-Release_1_8_20" ]; then
  wget -O - https://github.com/doxygen/doxygen/archive/Release_1_8_20.tar.gz | tar xz
  cd doxygen-Release_1_8_20
  mkdir build
  cd build
  cmake -Duse_libclang=ON  ..
  make VERBOSE=1 -j2
fi
export PATH=${DEPS_DIR}/doxygen-Release_1_8_20/build/bin:$PATH
doxygen --version

#Change back to where we left off.
cd $dir
