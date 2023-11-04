#!/bin/bash
#SBATCH --qos=regular
#SBATCH --time=0:30:00
#SBATCH --nodes=18
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=128
#SBATCH --constraint=cpu
#SBATCH --mem=400GB
#SBATCH --mail-type=END
#SBATCH --mail-user=yuxi.hong.research@gmail.com

export OMP_NUM_THREADS=128
srun /global/homes/y/yuxihong/graphclustering/CombBLAS/build/ReleaseTests/ReadWriteMtx \
/pscratch/sd/y/yuxihong/hipmcl_dataset/portal/Renamed_arch_vs_arch_30_50length.indexed.mtx \
/pscratch/sd/y/yuxihong/hipmcl_dataset/portal/arch_vs_arch_30_50length_propermm.mtx 0 

