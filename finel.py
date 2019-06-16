import numpy as np
import matplotlib.pyplot as plt
'''
Solve the differential equation:

x' + x = 0, x'(1) = 1

Usually differential equation solvers don't have the ability to 

x(t) = -exp(1-t)

'''

tfinal = 1.0
xdotfinal = 1.0

N = 100
A = np.zeros((N,N))
dx = tfinal/N

for i in range(1,N-1):
    A[i][i] = dx
    A[i][i-1] = -0.5
    A[i][i+1] = 0.5

A[0][0] = (dx-(3/2))
A[0][1] = 2
A[0][2] = -1/2
A[-1][-1] = 3/2
A[-1][-2] = -2
A[-1][-3] = 1/2
b = np.zeros(N)
b[-1] = xdotfinal*dx

soln = np.linalg.solve(A,b)
t = np.linspace(0,1,N)
real = -1*np.exp(1-t)
plt.plot(t,soln)
plt.plot(t,real)
plt.show()
#======================================================================
tfinal = 1.0
xdotfinal = 1.0

N = 100
A = np.zeros((N,N))
dx = tfinal/N

for i in range(1,N-1):
    A[i][i] = 1
    A[i][i-1] = -0.5/dx
    A[i][i+1] = 0.5/dx

A[0][0] = (dx-(3/2))
A[0][1] = 2
A[0][2] = -1/2
A[-1][-1] = 3/2
A[-1][-2] = -2
A[-1][-3] = 1/2
b = np.zeros(N)
b[-1] = xdotfinal*dx

soln = np.linalg.solve(A,b)
t = np.linspace(0,1,N)
real = -1*np.exp(1-t)
plt.plot(t,soln)
plt.plot(t,real)
plt.show()

rms = np.average(np.sqrt((real-soln)**2))
