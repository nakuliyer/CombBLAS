#!/bin/bash
#SBATCH --qos=debug
#SBATCH --time=0:05:00
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=128
#SBATCH --constraint=cpu
#SBATCH --mem=400GB
#SBATCH --mail-type=END
#SBATCH --mail-user=yuxi.hong.research@gmail.com

# grid 2x2 = 4 mpi process
# mpi process per node 2
# node 2
export OMP_NUM_THREADS=128
SUBMITROOT=/pscratch/sd/y/yuxihong/graphclustering/CombBLAS_submit
srun $SUBMITROOT/build/Applications/spgemm1d $VIRUS
# srun /global/homes/y/yuxihong/graphclustering/CombBLAS/build/Applications/spgemm1d $ARCH
# srun /global/homes/y/yuxihong/graphclustering/CombBLAS/build/Applications/spgemm1d $EUK
# srun /global/homes/y/yuxihong/graphclustering/CombBLAS/a.out
