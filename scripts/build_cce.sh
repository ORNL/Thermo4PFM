#!/bin/csh
module load PrgEnv-cray
module load craype-accel-amd-gfx90a
module load boost
module load cmake
module load cray-python
module load rocm

set ROOT = `pwd`

set BUILD_DIR = ${ROOT}/build
rm -rf ${BUILD_DIR}

mkdir -p ${BUILD_DIR}
cd ${BUILD_DIR}

# call cmake
cmake -DCMAKE_CXX_COMPILER=CC \
      -DCMAKE_BUILD_TYPE=Release \
      -DWITH_OPENMP_OFFLOAD=ON \
      -DCMAKE_CXX_FLAGS="-fopenmp" \
      -DMPIEXEC_EXECUTABLE="srun" \
      -DMPIEXEC_NUMPROCS_FLAG="-n" \
      -DMPIEXEC_NUMPROCS="1" \
      -DMPIEXEC_PREFLAGS="-c1;--gpus=1" \
      -DCMAKE_EXE_LINKER_FLAGS="-lm -fopenmp" \
      ..

make
