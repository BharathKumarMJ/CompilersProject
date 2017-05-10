import scipy
from scipy.io import mmread, mmwrite
from scipy import sparse
from scipy.sparse.linalg import LinearOperator
import numpy as np


N = 10000
sparsity = 0.2
denseMat = np.random.rand(N,N)
for i in range (N):
	for j in range (N):
		denseMat[i][j] = denseMat[i][j] * np.random.binomial(1,sparsity)
print "Made dense matrix"
S = sparse.csr_matrix(denseMat)
print "Generated csr matrix.  Writing to file..."
mmwrite("synth_dense_"+str(N)+".mtx", S)
print "File IO complete"
