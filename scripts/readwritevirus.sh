#!/bin/bash
#SBATCH --qos=debug
#SBATCH --time=0:05:00
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=128
#SBATCH --constraint=cpu

# grid 2x2 = 4 mpi process
# mpi process per node 2 
# node 2
module unload cray-mpich
module unload cray-libsci
module use /global/common/software/m3169/perlmutter/modulefiles
module load openmpi
mpicc /global/homes/y/yuxihong/graphclustering/CombBLAS/helloworld.c
# . /global/common/software/m4293/intelmkl/mkl/2023.2.0/env/vars.sh   
export OMP_NUM_THREADS=128
mpirun -np 4 /global/homes/y/yuxihong/graphclustering/CombBLAS/a.out
# srun --mpi=pmix -n 4 /global/homes/y/yuxihong/graphclustering/CombBLAS/build/ReleaseTests/ReadWriteMtx \
# /pscratch/sd/y/yuxihong/hipmcl_dataset/portal/Renamed_vir_vs_vir_30_50length.indexed.mtx \
# /pscratch/sd/y/yuxihong/hipmcl_dataset/portal/vir_vs_vir_30_50length_propermm.mtx 0 

