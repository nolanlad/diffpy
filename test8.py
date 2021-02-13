import numpy as np 
import matplotlib.pyplot as plt 
from utilities import *

'''
Solve for steady state solution for heat conduction in a
square of material with adiabatic walls and constant temperatures
at top and bottom
'''
N=20
x = np.linspace(0,1,N)
y = np.linspace(0,1,N)

xx,yy = np.meshgrid(x,y)
xx = xx.reshape(-1)
yy = yy.reshape(-1)



from Node2 import *

n = Nodes({"x":xx,"y":yy})

edgex1 = n.get_dim('x') == 1.0
edgex2 = n.get_dim('x') == 0.0
edgey1 = n.get_dim('y') == 1.0
edgey2 = n.get_dim('y') == 0.0
labels = np.zeros(n.lennodes)
labels[edgex1] = 1; labels[edgex2] = 2
labels[edgey1] = 3; labels[edgey2] = 4
labnam = {"mass":0,"edgex1":1,"edgex2":2,"edgey1":3,"edgey2":4}
n.add_labels(labels,labnam)
# n.plot('x','y')
# plt.show()

h = 3 # for some reason h < 3 produces unstable solutions

Axx = make_stiffness(n,2,h,'x')
Ayy = make_stiffness(n,2,h,'y')
Ax =  make_stiffness(n,1,h,'x')
AI =  make_stiffness(n,0,h,'x')
Ay =  make_stiffness(n,1,h,'y')

A = np.zeros(Axx.shape)

Axxyy = Axx + Ayy

mass = np.where(n.labels == n.labnames['mass'])[0]

for i in mass:
    A[i] = Axxyy[i]

isedgex = (n.labels == n.labnames['edgex1'])|(n.labels == n.labnames['edgex2'])


edgex = np.where(isedgex)[0]

for i in edgex:
    A[i] = AI[i]

isedgey = (n.labels == n.labnames['edgey1'])|(n.labels == n.labnames['edgey2'])
edgey = np.where(isedgey&(~isedgex))[0]

for i in edgey:
    A[i] = AI[i]

edgex1 = np.where((n.labels == n.labnames['edgey2']))[0]

b = np.zeros(n.lennodes)
# b = np.ones(n.lennodes)
for i in edgex1:
    b[i] = 100

edgex2 = np.where((n.labels == n.labnames['edgex2']))[0]

for i in edgex2:
    b[i] = 100

edgex1 = np.where((n.labels == n.labnames['edgex1']))[0]

# b = np.ones(n.lennodes)
for i in edgex1:
    b[i] = 400

edgex2 = np.where((n.labels == n.labnames['edgey1']))[0]

for i in edgex2:
    b[i] = 400

T = svdsolve(A,b)
Ts = T.reshape((N,N))
plt.imshow(Ts)
plt.colorbar()
plt.show()

import numpy.ma as ma

plt.imshow(ma.masked_equal(A,0))
plt.show()
print(np.linalg.det(A))