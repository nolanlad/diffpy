'''
Steady state temperature in a 1 meter aluminum bar with steady 
heat generation of 1W/m^3
'''

from diff import *
import matplotlib.pyplot as plt


#thermal properties of Aluminum
k = 205.
alpha = 9.7e-5
Q = 1.0

N = 25
x = np.sqrt(np.linspace(0,1,N))
nods = nodelinspace_1D(x)
A = np.zeros((N,N))

for i in range(1,N-1):
	ids,vals = nods[i].diff('x',2,0)
	A[i][ids] += alpha*vals

#constant temperature boundary condition
A[0][0] = 1
#A[-1][-1] = 1

#Adiabatic
ids,vals = nods[-1].diff('x',1,0)
A[-1][ids] = -k*vals

#boundary conditions
b = -Q*np.ones(N)/k
b[0] = 300 #kelvin
#b[-1] = 300 #kelvin
b[-1] = 0
T = np.linalg.solve(A,b)

print('max temp = ', np.amax(T))
plt.plot(x,T,'bo')
plt.show()

