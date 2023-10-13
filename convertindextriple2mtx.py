# %%
import numpy as np
from os.path import join, exists
import sys
if len(sys.argv) < 2:
    print("usage: python convertindextriple2mtx.py index.triple.file")
    sys.exit()
infile = sys.argv[1]
# infile = "/Users/hongy0a/Documents/dataset/HipMCL/vir_vs_vir_30_50length.indexed.triples"
outfile = infile + '.mtx'

# read distmapper.mtx 
IndexTripleFormat=infile
m = {}
index = 0
outstr = []
with open(IndexTripleFormat,'r') as f:
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
        outstr.append( "{}\t{}\t{}\n".format(m[n1],m[n2],val) )
# %%
# write index.triple format to mtx format.
wf = open(outfile,'w')
for line in outstr:
    wf.write(line)
wf.close()
# %%
