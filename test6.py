from diff import *
import matplotlib.pyplot as plt

noder = node_rect()
#thermal properties of Aluminum
k = 205.
alpha = 9.7e-5
Q = 1.0

noder2 = []
for line in noder:
    noder2.extend(line)

N = len(noder2)
A = np.zeros((N,N))
for i in range(N):
    ids,vals = noder2[i].diff('x',1,2)
    A[i][ids] += alpha*vals

x = np.array([n.get_val('x') for n in noder2])
y = np.array([n.get_val('y') for n in noder2])

z = (x**2)*(y**2)
dzdx = (2*x*y**2)

dzdxa = np.matmul(A,z)
print('f word')
plt.plot(dzdx,dzdxa,'ro')
plt.show()