import scipy
from scipy.io import mmread, mmwrite
from scipy import sparse
from scipy.sparse.linalg import LinearOperator
import numpy as np
# import matplotlib.pyplot as plt
import random

N = 131072
block = 4
test = 0

csrMat = sparse.lil_matrix((N,N))
filename = "synth_block_"+str(N)+"_"+str(block)+".mtx"
# file = open(filename, 'w')
# file.write("%%MatrixMarket matrix coordinate real general\n%\n")
# file.write(str(N) + " " + str(N) + " " + str(N*block) + "\n")
for i in range (N):
	if(i % block == 0):
		for ii in range (block):  
			for jj in range (block):
				# file.write(str(ii+i+1) + " " + str(jj+i+1) + " " + str(random.random()) + "\n")
				# file.write(str(ii+i+1) + " " + str((jj+i+ N/2)%N + 1) + " " + str(random.random()) + "\n")
				csrMat[ii+i, (jj+i+N/2)%N] = random.random()
				csrMat[ii+i, jj+i] = random.random()
mmwrite(filename, csrMat)
print "File IO complete"

# if test==1:
# 	# print denseMat
# 	# plt.imshow(denseMat)
# 	plt.spy(csrMat)	
# 	plt.show()
