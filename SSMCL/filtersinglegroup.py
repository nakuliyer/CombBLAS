# filter the group that only contains itself

import numpy as np
from os.path import join, exists
import sys
from os import remove
from collections import defaultdict

if len(sys.argv) < 2:
    print("missing args")
    sys.exit()

sg = 0
with open(sys.argv[1], 'r') as f:
    while 1:
        line = f.readline()
        if not line:
            break
        sl = line.strip('\n').split(' ')
        if len(sl) == 1:
            sg += 1
print("single group size ", sg)


totalnodes = []
with open(sys.argv[1], 'r') as f:
    while 1:
        line = f.readline()
        if not line:
            break
        sl = line.strip('\n').split(' ')
        isl = []
        for x in sl:
            try:
                isl.append(int(x))
            except:
                pass
        totalnodes.extend(isl)
print("totalnodes size ", len(set(totalnodes)))
print(sorted(totalnodes)[:10])
print(sorted(totalnodes)[-10:])


totalnodes = []
with open('vir_vs_vir_30_50length.indexed.triples.mtx', 'r') as f:
    line = f.readline()
    line = f.readline()
    while 1:
        line = f.readline()
        if not line:
            break
        sl = line.strip('\n').split('\t')[:2]
        isl = []
        for x in sl:
            try:
                isl.append(int(x))
            except:
                pass
        totalnodes.extend(isl)
print("totalnodes size ", len(set(totalnodes)))
print(sorted(totalnodes)[:10])
print(sorted(totalnodes)[-10:])