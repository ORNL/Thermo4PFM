#!/bin/csh -f
module load gcc/11.1.0
module load boost
module load cmake
module load cuda
module load python

set ROOT = `pwd`

set BUILD_DIR = ${ROOT}/build
mkdir -p ${BUILD_DIR}
cd ${BUILD_DIR}

# call cmake
cmake -DCMAKE_CXX_COMPILER=mpicxx \
      -DCMAKE_BUILD_TYPE=Release \
      -DMPIEXEC_EXECUTABLE="/sw/summit/xalt/1.2.1/bin/jsrun" \
      -DMPIEXEC_NUMPROCS_FLAG="-n" \
      -DMPIEXEC_NUMPROCS="1" \
      -DMPIEXEC_PREFLAGS="-a1;-c4;-bpacked:2;-g1" \
      -DWITH_OPENMP_OFFLOAD=ON \
      -DCMAKE_CXX_FLAGS="-fopenmp -foffload=-lm -foffload=nvptx-none -foffload="-O3" -fno-stack-protector" \
      -DCMAKE_EXE_LINKER_FLAGS="-lm -fopenmp" \
      ..

make -j 4
