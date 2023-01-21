#!/bin/csh -f
#SBATCH -p batch
#SBATCH -A MAT190_crusher
#SBATCH -J TernaryThermo4PFMcrusher
#SBATCH -o TernaryThermo4PFMcrusher.o%J
#SBATCH -t 00:15:00
#SBATCH -N 1
module use /sw/crusher/ums/compilers/modulefiles
module load llvm

module load boost
module load craype-accel-amd-gfx90a
module load rocm

cd build

ln -fs ../thermodynamic_data/calphadMoNbTa.json

setenv OMP_NUM_THREADS 1
setenv LIBOMPTARGET_INFO 31

set exec = ./drivers/loopCALPHADConcSolverTernary
#srun -n1 -N1 -c1 --gpus=1 $exec 1000
#srun -n1 -N1 -c1 --gpus=1 $exec 10000
srun -n1 -N1 -c1 --gpus=1 $exec 100000
#srun -n1 -N1 -c1 --gpus=1 $exec 1000000
#srun -n1 -N1 -c1 --gpus=1 $exec 10000000

