#/bin/bash

ROOT=`pwd`

BUILD_DIR=${ROOT}/build
mkdir -p ${BUILD_DIR}
cd ${BUILD_DIR}

# call cmake
cmake -DCMAKE_CXX_COMPILER=mpiCC \
      -DWITH_CLANG_FORMAT=ON \
      -DCMAKE_INSTALL_PREFIX=${HOME}/Thermo4PFM \
      ..

make -j 4
