import scipy
from scipy.io import mmread, mmwrite
from scipy import sparse
from scipy.sparse.linalg import LinearOperator
import numpy as np
import matplotlib.pyplot as plt
import random

N = 1048576
block = 16
# cooMat = sparse.coo_matrix((N,N))
filename = "synth_block_"+str(N)+"_"+str(block)+".mtx"
file = open(filename, 'w')
file.write("%%MatrixMarket matrix coordinate real general\n%\n")
file.write(str(N) + " " + str(N) + " " + str(N*block) + "\n")
for i in range (N):
	if(i % block == 0):
		for ii in range (block):
			for jj in range (block):
				file.write(str(ii+i+1) + " " + str(jj+i+1) + " " + str(random.random()) + "\n")
print "File IO complete"

# print denseMat
# plt.imshow(denseMat)
# plt.show()
