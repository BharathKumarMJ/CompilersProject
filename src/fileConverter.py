import scipy
from scipy.io import mmread, mmwrite
from scipy.sparse.linalg import LinearOperator
import numpy as np
import sys

N = int(sys.argv[2])	

filename = sys.argv[1]
S = mmread(filename)

B = S.tobsr((N,N))

print "Done converting to bsr"

data = B.data
indices = B.indices
indptr = B.indptr
bm, bn = np.asarray(B.get_shape())/N
r, c = N, N
nnzb = data.shape[0]

outFile = open(filename.split(".mtx")[0]+"_"+str(N), 'w')

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
outFile.write("% Data: \n")
for i in range(data.size/(N*N)):
	for j in range(N):
		for k in range(N):
			outFile.write(str( data[i][j][k] ) + "\n")
	#outFile.write("\n")
	# print str( (data.flatten())[i] )

outFile.write("END \n")
outFile.close()
# print B.matmat(B)