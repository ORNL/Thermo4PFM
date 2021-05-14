#/bin/bash
module load PE-gnu
module load boost
module load cmake
module load python

ROOT=`pwd`

BUILD_DIR=${ROOT}/build
mkdir -p ${BUILD_DIR}
cd ${BUILD_DIR}

# call cmake
cmake -DCMAKE_CXX_COMPILER=mpiCC \
      -DWITH_CLANG_FORMAT=ON \
      -DCMAKE_PREFIX_PATH=${HOME}/bin \
      ..

make -j 4
