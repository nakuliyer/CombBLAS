# %%
######################################
# Import Module and Function Definition
######################################
import sys
import numpy as np
from os.path import join, exists
import matplotlib.pyplot as plt
from scipy.sparse import csr_array 

def loadmtx(filename):
    ridx = []
    cidx = []
    values = []
    with open(filename,'r') as f:
        f.readline()
        m,n,nnzval = f.readline().strip('\n').split(' ')
        m,n,nnzval = int(m),int(n),int(nnzval)
        print(m,n,nnzval)
        nnzcnt = 0
        while 1:
            line = f.readline()
            if not line:
                break
            try:
                n1,n2,val = line.strip('\n').split('\t')
            except:
                print("error ", line)
                sys.exit()
            ridx.append(int(n1))
            cidx.append(int(n2))
            values.append(float(val))
            nnzcnt += 1
        assert(nnzval == nnzcnt)
    return csr_array((values,(ridx,cidx)),shape=(m,n),dtype=np.float64)

def checkequal(spA,spB):
    #TODO: improve it
    return np.allclose(spA.indptr,spB.indptr) and \
    np.allclose(spA.indices,spB.indices) and np.allclose(spA.data, spB.data)


def partial_sort_with_indices(arr, start, end):
    # Enumerate the array elements with their indices
    indexed_arr = list(enumerate(arr[start:end], start=start))
    # Sort the enumerated array based on element value
    sorted_arr = sorted(indexed_arr, key=lambda x: x[1])
    # Extract and return the sorted indices
    sorted_indices = [index for index, _ in sorted_arr]
    return sorted_indices

def sortedscipysparsearray(sparray):
    for i in range(len(sparray.indptr)-1):
        if sparray.indptr[i] == sparray.indptr[i+1]:
            continue
        startindices = sparray.indptr[i]
        endindices = sparray.indptr[i+1]
        sorted_indices = partial_sort_with_indices(sparray.indices, startindices, endindices)
        sparray.indices[startindices:endindices] = sparray.indices[sorted_indices]
        sparray.data[startindices:endindices] = sparray.data[sorted_indices]
    return sparray

# %%
IndexTripleFormat="/Users/hongy0a/Documents/dataset/HipMCL/vir_vs_vir_30_50length.indexed.triples"
# IndexTripleFormat = "/Users/hongy0a/Documents/Research/graphclustering/CombBLAS/Test1k"
sparray = loadmtx(IndexTripleFormat + ".mtx")
print(sparray.shape)
# sparray_read = loadmtx(IndexTripleFormat + "_reading.mtx")
# print(checkequal(sparray_read,sparray))

# %%
# plt.figure()
# plt.imshow(sparray[:1000,:1000].todense())
# plt.show()
# %%
# remove diag 
BS = 54929
def diagop(sparray,keepdiag):
    m,n = sparray.shape
    ridx = []
    cidx = []
    values = []
    for ri in range(len(sparray.indptr)-1):
        startidx = sparray.indptr[ri]
        endidx = sparray.indptr[ri+1]
        for idx in range(startidx,endidx):
            cj = sparray.indices[idx]
            val = sparray.data[idx]
            op = False
            if keepdiag:
                op = int(ri / BS) == int(cj / BS)
            else:
                op = int(ri / BS) != int(cj / BS)
            if op:
                ridx.append(ri)
                cidx.append(cj)
                values.append(val)
    return csr_array((values,(ridx,cidx)),shape=(m,n),dtype=np.float64)
# spoffdiag = diagop(sparray,keepdiag=False)
# spoffdiag_cpp = loadmtx(IndexTripleFormat + "_offdiag.mtx")
# print(checkequal(spoffdiag,spoffdiag_cpp))
# spdiag = diagop(sparray,keepdiag=True)
# spdiag_cpp = loadmtx(IndexTripleFormat + "_diag.mtx")
# print(checkequal(spdiag,spdiag_cpp))

# %%
spA = loadmtx(IndexTripleFormat + ".mtx")
spB = loadmtx(IndexTripleFormat + ".mtx")
spA = diagop(spA,keepdiag=True)
spB = diagop(spB,keepdiag=True)
spC = spA @ spB

# %%
spC2d = loadmtx(IndexTripleFormat+"_C.mtx")
spCpy_sorted = sortedscipysparsearray(spC)
spC2d_sorted = sortedscipysparsearray(spC2d)
print("check 2D gemm",checkequal(spCpy_sorted,spC2d_sorted))
# %%
plt.figure()
plt.imshow(spC2d[:1000,:1000].todense())
plt.show()
# %%
spC1d = loadmtx(IndexTripleFormat+"_From1D.mtx")
print(spC1d.shape)
spC1d_sorted = sortedscipysparsearray(spC1d)
print("check 1d gemm",checkequal(spCpy_sorted,spC1d_sorted))

# %%
plt.figure()
plt.imshow(spC1d[:1000,:1000].todense())
plt.show()
# %%
def checkparts(spA, spB, N):
    print(checkequal(sortedscipysparsearray(spA[:N,:N]),sortedscipysparsearray(spB[:N,:N])))
# %%
checkparts(spC2d,spC1d,BS+1)
# %%
print(spC1d.shape)
print(spC2d.shape)
# plt.figure()
# plt.imshow(spC1d[,BS:(BS+1000)].todense())
# plt.show()
# %%
