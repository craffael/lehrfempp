# This script installs an up-2-date version of cmake on the github actions build server.

# Exit immediately from this script upon error
set -e

# save current directory
dir=$(pwd)

mkdir -p ${DEPS_DIR} && cd ${DEPS_DIR}

if [ ! -d "doxygen-1.9.1" ]; then
  wget -O - https://sourceforge.net/projects/doxygen/files/rel-1.9.1/doxygen-1.9.1.src.tar.gz/download | tar xz
  cd doxygen-1.9.1
  mkdir build
  cd build
  export CXX=g++-9
  cmake -Duse_libclang=ON  ..
  make -j2
fi
export PATH=${DEPS_DIR}/doxygen-1.9.1/build/bin:$PATH
doxygen --version

#Change back to where we left off.
cd $dir
