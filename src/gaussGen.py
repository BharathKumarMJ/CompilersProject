import numpy as np
import matplotlib.pyplot as plt
from random import randint

numSamples = 10
numIters = 1
mean = [0, 0]
X = np.zeros((numSamples, numSamples))
cov = [[(numSamples/numIters), 0], [0, (numSamples/numIters)]]  # diagonal covariance
for k in range(numIters):
	mean = [(numSamples/numIters)*k, (numSamples/numIters)*k]
	# [randint(-20,20)*np.sum(np.random.rand(1)),randint(-20,20)*np.sum(np.random.rand(1))]
	C = np.random.multivariate_normal(mean, cov, numSamples/numIters).T
	# plt.plot(C[0], C[1], 'x')
	for i in range(C[0].shape[0]):
		for j in range(C[1].shape[0]):
			II = min(max(0, np.rint(C[0][i]) ), X.shape[0]-1)
			JJ = min(max(0, np.rint(C[1][j]) ), X.shape[1]-1)
			X[II][JJ] += 100*np.random.random(1)
			print II, JJ, X[II][JJ]

print X
plt.imshow(X)
# plt.colorbar(orientation='vertical')

# exit()
# # print X.shape
		# print II, JJ, X[II][JJ]
plt.axis([0, 4096, 0, 4096])
plt.show()