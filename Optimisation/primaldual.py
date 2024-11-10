# %%
import numpy as np
import scipy


n = 2
constraints = 1
A = scipy.linalg.hilbert(n)
B = np.ones(n)
c = np.array(1)
b = np.dot(A, np.ones(n))

Ainv = np.linalg.inv(A)
p = (np.dot(np.dot(B, Ainv), b) - c)/(np.dot(B, np.dot(Ainv,B.T)))
x = np.dot(Ainv, (b - np.dot(B.T, p)))
print(x,p)
print(x.shape, p.shape)
# BA−1Btp = BA−1b − c,
