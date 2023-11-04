# %%
import numpy as np
from os.path import join, exists
import sys
from os import remove

if len(sys.argv) < 2:
    print("usage: python convertindextriple2mtx.py index.triple.file")
    sys.exit()
infile = sys.argv[1]
outfile = infile + '.mtx'
if exists(outfile):
    remove(outfile)

# read index.triple format file
m = {}
matrixmap = {}
index = 0
outstr = []
nnz = 0
onebased = True #!!! ALWAYS BE TRUE
if onebased:
    index += 1
maxr,maxc = 0,0
with open(infile,'r') as f:
    while 1:
        line = f.readline()
        if not line:
            break
        n1,n2,val = line.strip('\n').split('\t')
        if n1 not in m:
            m[n1] = index
            index += 1
        if n2 not in m:
            m[n2] = index
            index += 1
        n1 = m[n1]
        n2 = m[n2]
        maxr = max(maxr,n1)
        maxc = max(maxc,n2)
        nnz += 1
        matrixmap['{}_{}'.format(n1,n2)] = val
        outstr.append( "{}\t{}\t{}\n".format(n1,n2,val) )
assert(maxr == maxc)
if not onebased:
    maxr += 1
    maxc += 1
# %%
# write index.triple format to mtx format.
with open(outfile, 'w') as wf:
    wf.write("%%MatrixMarket matrix coordinate real general\n")
    wf.write("{} {} {}\n".format(maxr,maxc,nnz))
    for line in outstr:
        wf.write(line)


#TODO: force input matrix to be symmetric, need to verify the correctness
for i in range(1,maxr+1): # one based
    for j in range(1,maxc+1): 
        indexstr = '{}_{}'.format(i,j)
        if indexstr in matrixmap:
            transindexstr = '{}_{}'.format(j,i)
            if transindexstr in matrixmap:
                maxval = max(matrixmap[indexstr],matrixmap[transindexstr])
            else:
                matrixmap[transindexstr] = matrixmap[indexstr]
outstr = []
for i in range(1,maxr+1): # one based
    for j in range(1,maxc+1): 
        indexstr = '{}_{}'.format(i,j)
        if indexstr in matrixmap:
            outstr.append( "{}\t{}\t{}\n".format(i,j,matrixmap[indexstr]) )
if exists(outfile):
    remove(outfile)
# write index.triple format to mtx format.
with open(outfile, 'w') as wf:
    wf.write("%%MatrixMarket matrix coordinate real general\n")
    wf.write("{} {} {}\n".format(maxr,maxc,nnz))
    for line in outstr:
        wf.write(line)
