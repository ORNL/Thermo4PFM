#!/bin/csh -f
module use /sw/summit/modulefiles/ums/stf010/Core/
module load llvm/14

module load boost
module load cmake
module load python
module load cuda

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
      -DCMAKE_CXX_FLAGS="-fopenmp -fopenmp-targets=nvptx64 -Rpass=openmp-opt -Rpass-missed=openmp-opt" \
      -DWITH_OPENMP_OFFLOAD=ON \
      -DMPIEXEC_EXECUTABLE="srun" \
      -DMPIEXEC_NUMPROCS_FLAG="-n" \
      -DMPIEXEC_NUMPROCS="1" \
      -DCMAKE_EXE_LINKER_FLAGS="-lm -fopenmp" \
      ..

make
