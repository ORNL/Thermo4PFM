#!/bin/bash
#BSUB -q debug
#BSUB -P MAT190
#BSUB -J LLVMTernaryThermo4PFM
#BSUB -o LLVMTernaryThermo4PFM.o%J
#BSUB -W 1:00
#BSUB -nnodes 1
#BSUB -env "all"
module use /sw/summit/modulefiles/ums/stf010/Core/
module load llvm/15

module load boost
module load cmake
module load python
module load cuda

cd build

ln -fs ../thermodynamic_data/calphadMoNbTa.json

export OMP_NUM_THREADS=1
##export LIBOMPTARGET_INFO=31

###jsrun -n6 -a1 -c7 -g1 -d packed -b packed:7 nvprof --print-gpu-trace ./drivers/loopCALPHADConcSolverTernary 1000000
###jsrun -n6 -a1 -c7 -g1 -d packed -b packed:7 ./drivers/loopCALPHADConcSolverTernary 1000000

##exec=../bin/loopCALPHADConcSolverTernary_nolto
exec=./drivers/loopCALPHADConcSolverTernary
nreps=2
jsrun -n1 -a1 -c7 -g1 -d packed -b packed:7 nvprof $exec 1000 $nreps
jsrun -n1 -a1 -c7 -g1 -d packed -b packed:7 nvprof $exec 10000 $nreps
jsrun -n1 -a1 -c7 -g1 -d packed -b packed:7 nvprof $exec 100000 $nreps
jsrun -n1 -a1 -c7 -g1 -d packed -b packed:7 nvprof $exec 1000000 $nreps
jsrun -n1 -a1 -c7 -g1 -d packed -b packed:7 nvprof $exec 10000000 $nreps


#export OMP_NUM_THREADS=4
#jsrun -n1 -a1 -c7 -g1 -d packed -b packed:7 ./drivers/loopCALPHADConcSolverTernary $n

#export OMP_NUM_THREADS=8
#jsrun -n1 -a1 -c8 -g1 -d packed -b packed:8 ./drivers/loopCALPHADConcSolverTernary $n

