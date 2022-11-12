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
      -DCMAKE_BUILD_TYPE=Release \
      -DWITH_CLANG_FORMAT=ON \
      -DWITH_CONVERGENCE_HISTORY=ON \
      -DCMAKE_INSTALL_PREFIX=${HOME}/Thermo4PFM \
      -DCMAKE_PREFIX_PATH=${HOME}/bin \
      ..

make -j 4
