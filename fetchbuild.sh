#!/bin/bash
#
ROOTF=/pscratch/sd/y/yuxihong/graphclustering/CombBLAS_submit
rsync -avrz --exclude "build" --exclude ".git" --exclude ".cache" /pscratch/sd/y/yuxihong/graphclustering/CombBLAS/ /pscratch/sd/y/yuxihong/graphclustering/CombBLAS_submit
cd $ROOTF/build
make -j
cd $ROOTF

