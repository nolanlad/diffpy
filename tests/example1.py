'''
Steady state temperature in a 1 meter aluminum bar with steady 
heat generation of 1W/m^3
'''

from diff import *
import matplotlib.pyplot as plt
import scipy.sparse
import scipy


#thermal properties of Aluminum
k = 205.
alpha = 9.7e-5
Q = 1.0

N = 10
accuracy = 4
x = np.sqrt(np.linspace(0,1,N))
nods = nodelinspace_1D(x)
A = np.zeros((N,N))


for i in range(1,N-1):
	ids,vals = nods[i].diff('x',2,accuracy)
	A[i][ids] += alpha*vals

#constant temperature boundary condition
A[0][0] = 1
#A[-1][-1] = 1

#Adiabatic
ids,vals = nods[-1].diff('x',1,accuracy)
A[-1] = 0
A[-1][ids] = -k*vals

#boundary conditions
b = -Q*np.ones(N)/k
b[0] = 300 #kelvin
#b[-1] = 300 #kelvin
b[-1] = 0
T = np.linalg.solve(A,b)

print('max temp = ', np.amax(T))
plt.plot(x,T,'b-')
#plt.show()

ids,vals = nods[0].diff('x',1,accuracy)
A[0] = 0
A[0][ids] =  -alpha*vals
ids,vals = nods[-1].diff('x',1,accuracy)
A[-1] = 0
A[-1][ids] =  -alpha*vals

A = scipy.sparse.csr_matrix(A)

dt = 0.1
T2 = T
for i in range(10000):
    #dT = dt*np.matmul(A,T2)
    dT = dt*(A@T2)
    T2 = T2 + dT
    T2[0] = 300
    T2[-1] = T2[-2]

plt.plot(x,T2,'r-')

for i in range(30000):
    #dT = dt*np.matmul(A,T2)
    dT = dt*(A@T2)
    T2 = T2 + dT
    T2[0] = 300
    T2[-1] = T2[-2]

print(np.amax(T2))

plt.plot(x,T2,'g-')