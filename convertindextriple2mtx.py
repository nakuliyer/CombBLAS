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
index = 0
outstr = []
nnz = 0
onebased = False
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
        outstr.append( "{}\t{}\t{}\n".format(n1,n2,val) )
# we are not onebased
maxr += 1
maxc += 1
# %%
# write index.triple format to mtx format.
with open(outfile, 'w') as wf:
    wf.write("%%MatrixMarket matrix coordinate real general\n")
    wf.write("{} {} {}\n".format(maxr,maxc,nnz))
    for line in outstr:
        wf.write(line)
