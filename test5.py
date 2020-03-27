import numpy as np 
import matplotlib.pyplot as plt 
from utilities import *


x = np.linspace(0,1,10)
y = np.linspace(0,1,10)

xx,yy = np.meshgrid(x,y)
xx = xx.reshape(-1)
yy = yy.reshape(-1)

x1 = xx[25]
y1 = yy[25]


from Node2 import *

n = Nodes({"x":xx,"y":yy})

i = 25
'''
M = n.lennodes
A = np.zeros((M,M))
Ay = np.zeros((M,M))
for i in range(1,n.lennodes):
    #print(i)
    # ids, soln = diff_x(i,n,1,1)
    ids, soln = diff_gen(i,n,1,1,'x')
    A[i][ids] = soln
    # ids, soln = diff_y(i,n,1,1)
    ids, soln = diff_gen(i,n,1,1,'y')
    Ay[i][ids] = soln
'''
A = make_stiffness(n,1,1,'x')
Ay = make_stiffness(n,1,1,'y')

# plt.imshow(Ay)
# plt.show()

z = (n.get_dim('x')**2)*n.get_dim('y')
dzdy = (n.get_dim('x')**2)
dza = np.matmul(A,z)
ddza = np.matmul(A,dza)
dz = 2*(n.get_dim('x')*n.get_dim('y'))
ddz =  2*n.get_dim('y')

# d = np.dot(z[ids],soln)

# plt.plot(ddza,ddz,'bo')
# plt.show()

dzdya = np.matmul(Ay,z)
# plt.plot(dzdya,dzdy,'bo')
# plt.show()

test1 = np.allclose(dzdya,dzdy)
test2 =  np.allclose(ddza,ddz)
test3 = np.allclose(dza,dz)

if test1 and test2 and test3:
    print('TEST 5 passed')


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

A = np.zeros(Axx.shape)

Axxyy = Axx + Ayy

mass = np.where(n.labels == n.labnames['mass'])[0]
print(mass)

for i in mass:
    A[i] = Axxyy[i]

edgex = np.where((n.labels == n.labnames['edgex1'])|(n.labels == n.labnames['edgex2']))[0]
print(edgex)

for i in edgex:
    A[i] = Ax[i]

edgey = np.where((n.labels == n.labnames['edgey1'])|(n.labels == n.labnames['edgey2']))[0]
print(edgey)

for i in edgey:
    A[i] = AI[i]

edgey1 = np.where((n.labels == n.labnames['edgey1']))[0]

b = np.zeros(n.lennodes)
for i in edgey1:
    b[i] = 300

edgey2 = np.where((n.labels == n.labnames['edgey2']))[0]

for i in edgey2:
    b[i] = 400

T = svdsolve(A,b)
Ts = T.reshape((10,10))
plt.imshow(Ts)
plt.colorbar()
plt.show()
plt.imshow(A)
plt.show()
'''
plt.imshow(Axxyy)
plt.show()

plt.imshow(Ax)
plt.show()

plt.imshow(AI)
plt.show()

plt.imshow(A)
plt.colorbar()
plt.show()
'''


