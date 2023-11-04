#!/bin/bash
#SBATCH --qos=regular
#SBATCH --time=0:30:00
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=128
#SBATCH --constraint=cpu
#SBATCH --mem=400GB

# grid 4x4 = 16 mpi process
# mpi process per node 2 
# node 8

export OMP_NUM_THREADS=128
SUBMITROOT=/pscratch/sd/y/yuxihong/graphclustering/CombBLAS_submit
srun $SUBMITROOT/build/Applications/spgemm1d $VIRUS
srun $SUBMITROOT/build/Applications/spgemm1d $ARCH
srun $SUBMITROOT/build/Applications/spgemm1d $EUK

