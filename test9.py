import numpy as np 
import matplotlib.pyplot as plt 
from utilities import *

'''
Solve for steady state solution for heat conduction in a
square of material with adiabatic walls and constant temperatures
at top and bottom
'''

N = 5

# x = np.linspace(0,1,N)
# y = np.linspace(0,1,N)
x = np.random.rand(N)
y = np.random.rand(N)
xx,yy = np.meshgrid(x,y)
xx = xx.reshape(-1)
yy = yy.reshape(-1)

xx = np.random.random(N**2)
yy = np.random.random(N**2)

from scipy.spatial import Delaunay
points = np.zeros((N*N,2))
points[:,0] = xx
points[:,1] = yy
tri = Delaunay(points)

print(points.shape)

plt.triplot(points[:,0], points[:,1], tri.simplices.copy())
plt.plot(points[:,0], points[:,1], 'o')
plt.show()