#!/bin/bash
#
ROOTF=/pscratch/sd/y/yuxihong/graphclustering/CombBLAS_submit
cd $ROOTF/build
rm -rf *
cmake ..
make -j
cd $ROOTF

