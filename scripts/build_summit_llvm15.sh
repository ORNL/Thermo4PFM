#!/bin/csh -f
module use /sw/summit/modulefiles/ums/stf010/Core/
module load llvm/15

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
      -DCMAKE_BUILD_TYPE=Release \
      -DWITH_OPENMP_OFFLOAD=ON \
      -DCMAKE_CXX_FLAGS="-fopenmp -foffload-lto -fopenmp-new-driver -fopenmp-cuda-mode -fopenmp-targets=nvptx64 -Rpass=openmp-opt -Rpass-missed=openmp-opt" \
      -DWITH_OPENMP_OFFLOAD=ON \
      -DMPIEXEC_EXECUTABLE="/sw/summit/xalt/1.2.1/bin/jsrun" \
      -DMPIEXEC_NUMPROCS_FLAG="-n" \
      -DMPIEXEC_NUMPROCS="1" \
      -DMPIEXEC_PREFLAGS="-a1;-c4;-bpacked:2;-g1" \
      -DCMAKE_EXE_LINKER_FLAGS="-lm -fopenmp" \
      ..

make
