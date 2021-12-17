#!/bin/csh -f
#setenv LLVM_ROOT /sw/spock/ums/eiw/llvm/14.0.0-20210819
setenv LLVM_ROOT /sw/spock/ums/eiw/llvm/14.0.0-20211207
setenv LD_LIBRARY_PATH ${LLVM_ROOT}/lib:/opt/rocm-4.2.0/lib:/opt/rocm-4.2.0/lib64:$LD_LIBRARY_PATH
setenv PATH ${LLVM_ROOT}/bin:$PATH

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
      -DCMAKE_BUILD_TYPE=Release \
      -DWITH_OPENMP_OFFLOAD=ON \
      -DCMAKE_CXX_FLAGS="-O2 -fopenmp -fopenmp-targets=amdgcn-amd-amdhsa -std=c++11 -D__NO_MATH_INLINES -U__SSE2_MATH__ -U__SSE_MATH__" \
      -DWITH_OPENMP_OFFLOAD=ON \
      -DMPIEXEC_EXECUTABLE="srun" \
      -DMPIEXEC_NUMPROCS_FLAG="-n" \
      -DCMAKE_EXE_LINKER_FLAGS="-lm -fopenmp" \
      ..

make
