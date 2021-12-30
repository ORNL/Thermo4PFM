#!/bin/csh -f
module load rocm-compiler/4.3.0
module load rocm/4.3.0
module load craype-accel-amd-gfx908

module load boost
module load cmake
module load cray-python

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
      -DCMAKE_CXX_FLAGS="-O3 -gline-tables-only -fopenmp -fopenmp-targets=amdgcn-amd-amdhsa -Xopenmp-target=amdgcn-amd-amdhsa -march=gfx908" \
      -DWITH_OPENMP_OFFLOAD=ON \
      -DMPIEXEC_EXECUTABLE="srun" \
      -DMPIEXEC_NUMPROCS_FLAG="-n" \
      -DMPIEXEC_NUMPROCS="1" \
      -DCMAKE_EXE_LINKER_FLAGS="-lm -fopenmp -fopenmp-targets=amdgcn-amd-amdhsa -Xopenmp-target=amdgcn-amd-amdhsa -march=gfx908" \
      ..

make
