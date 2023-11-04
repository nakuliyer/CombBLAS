import numpy as np
from os.path import join, exists
import sys
from os import remove
from collections import defaultdict

if len(sys.argv) < 2:
    print("usage: python convert2metisinput.py index.triple.file")
    sys.exit()

infile = sys.argv[1] + '.mtx'
outfile = sys.argv[1] + '.graph'
if exists(outfile):
    remove(outfile)
onebased = True
with open(infile,'r') as infilehandler:
    infilehandler.readline()
    m,n,nnz = infilehandler.readline().strip('\n').split(' ')
    maxnnz = 0.0
    nnzcount = 0
    while 1:
        line = infilehandler.readline()
        if not line:
            break
        n1,n2,val = line.strip('\n').split('\t')
        val = float(val)
        assert(val > 0.0)
        maxnnz = max(maxnnz, val)
        nnzcount +=1
    print("max value ", maxnnz, "nnz count ", nnzcount)

with open(outfile, 'w') as outfilehandler:
    with open(infile,'r') as infilehandler:
        infilehandler.readline()
        m,n,nnz = infilehandler.readline().strip('\n').split(' ')
        bigmap = defaultdict(str)
        nnzcount = 0
        while 1:
            line = infilehandler.readline()
            if not line:
                break
            n1, n2, val = line.strip('\n').split('\t')
            nnzcount += 1
            if not onebased:
                n1 = int(n1)
                n2 = int(n2)
                n1 += 1
                n2 += 1
                n1 = str(n1)
                n2 = str(n2) 
            val = float(val) * 10000
            val = str(int(val))
            # n1,n2,val = int(n1),int(n2),float(val)
            bigmap[n1] += n2 + " " + val + " "
            bigmap[n2] += n1 + " " + val + " "
    outfilehandler.write("{} {} 001\n".format(m,nnz))
    for i in range(1,int(m)+1):
        if str(i) in bigmap.keys():
            outfilehandler.write(" {}\n".format(bigmap[str(i)]))
        else:
            outfilehandler.write(" \n")

    assert('219715' in bigmap.keys())
    print("nnzcount", nnzcount)