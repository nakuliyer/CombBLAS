# %%
import numpy as np
from os.path import join, exists
import sys
from os import remove

if len(sys.argv) < 2:
    print("usage: python convertindextriple2mtx.py mtx.file")
    sys.exit()
infile = sys.argv[1]

# read index.triple format file
matrixmap = {}
outstr = []
maxr = maxc = 0
with open(infile,'r') as f:
    line = f.readline()
    maxr,maxc,nnz = f.readline().strip('\n').split(' ')
    maxr = int(maxr)
    maxc = int(maxc)
    while 1:
        line = f.readline()
        if not line:
            break
        n1,n2,val = line.strip('\n').split('\t')
        matrixmap['{}_{}'.format(n1,n2)] = val
assert(maxr == maxc)


#TODO: force input matrix to be symmetric, need to verify the correctness
for i in range(1,maxr+1): # one based
    for j in range(1,maxc+1): 
        indexstr = '{}_{}'.format(i,j)
        if indexstr in matrixmap:
            transindexstr = '{}_{}'.format(j,i)
            if transindexstr in matrixmap:
                assert(matrixmap[indexstr] == matrixmap[transindexstr])
            else:
                print("matrix not symm detected.")
                exit()
print("matrix is symm!")