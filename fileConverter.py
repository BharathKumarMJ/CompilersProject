import scipy
from scipy.io import mmread, mmwrite
from scipy.sparse.linalg import LinearOperator
import numpy as np

N = 25

S = mmread("tmt_unsym/tmt_unsym.mtx")

B = S.tobsr((N,N))
# print B.todense()
# print dir(B)

# print B.data
data = B.data
indices = B.indices
indptr = B.indptr
bm, bn = np.asarray(B.get_shape())/N
r, c = N, N
nnzb = data.shape[0]

outFile = open("tmt_unsym25", 'w')

# bn
outFile.write("% bn: \n")
outFile.write(str(bn))
outFile.write("\n")

# bm
outFile.write("% bm: \n")
outFile.write(str(bm))
outFile.write("\n")

# r
outFile.write("% R: \n")
outFile.write(str(r))
outFile.write("\n")

# c
outFile.write("% C: \n")
outFile.write(str(c))
outFile.write("\n")

# nnzb
outFile.write("% nnzb: \n")
outFile.write(str(nnzb))
outFile.write("\n")

# numIndices
outFile.write("% numIndices: \n")
outFile.write(str(indices.size))
outFile.write("\n")

# Indices
outFile.write("% Indices: \n")
for i in range(indices.size):
	outFile.write(str(indices[i]))
	outFile.write("\n")

# numPtrs
outFile.write("% numPtrs: \n")
outFile.write(str(indptr.size))
outFile.write("\n")

# Ptrs
outFile.write("% IndPtr: \n")
for i in range(indptr.size):
	outFile.write(str(indptr[i]))
	outFile.write("\n")

# Data size
outFile.write("% DataSize: \n")
outFile.write(str(data.size))
outFile.write("\n")

# Data
# print data
outFile.write("% Data: \n")
for i in range(data.size):
	outFile.write(str( (data.flatten())[i] ))
	outFile.write("\n")
	# print str( (data.flatten())[i] )

outFile.write("END \n")
outFile.close()
# print B.matmat(B)