#!/bin/csh -f
module use /sw/crusher/ums/compilers/modulefiles
module load llvm

module load boost
module load cmake
module load cray-python
module load craype-accel-amd-gfx90a
module load rocm

which cc
which clang
which clang++

set ROOT = `pwd`

set BUILD_DIR = ${ROOT}/build
mkdir -p ${BUILD_DIR}
cd ${BUILD_DIR}

# call cmake
cmake -DCMAKE_CXX_COMPILER=clang++ \
      -DCMAKE_BUILD_TYPE=Debug \
      -DWITH_OPENMP_OFFLOAD=ON \
      -DCMAKE_CXX_FLAGS="-fopenmp -foffload-lto -fopenmp -target x86_64-pc-linux-gnu -fopenmp-targets=amdgcn-amd-amdhsa -Xopenmp-target=amdgcn-amd-amdhsa -march=gfx90a" \
      -DWITH_OPENMP_OFFLOAD=ON \
      -DCMAKE_EXE_LINKER_FLAGS="-lm -fopenmp" \
      ..

make
