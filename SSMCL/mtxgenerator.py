import numpy as np
from os.path import join, exists
import sys
if len(sys.argv) < 5:
    # M is row, N is col, d is average nnz per row, filename is output name
    print("usage: python mtxgenerator.py M N d filename")
M = int(sys.argv[1])
N = int(sys.argv[2])
d = int(sys.argv[3])
filename = sys.argv[4]

# write index.triple format to mtx format.
with open(filename,'w') as f:
    f.write("%%MatrixMarket matrix coordinate real general\n")
    f.write("{} {} {}\n".format(M,N,d*M))
    # generate spmatrix
    for row in range(M):
        randcidx = np.random.choice(M, size=d, replace=False)
        randval = np.random.rand(d)
        for cidx,val in zip(randcidx,randval):
            f.write("{}\t{}\t{:.7f}\n".format(row,cidx,val))


