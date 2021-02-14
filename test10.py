import numpy as np 
import matplotlib.pyplot as plt 
from utilities import *

'''
Solve for steady state solution for heat conduction in a
square of material with adiabatic walls and constant temperatures
at top and bottom
'''
N=25
x = np.linspace(0,1,N)
y = np.linspace(0,1,N)

xx,yy = np.meshgrid(x,y)
xx = xx.reshape(-1)
yy = yy.reshape(-1)



from Node2 import *

n = Nodes({"x":xx,"y":yy})

edgex = (n.get_dim('x') == 1.0)|(n.get_dim('x') == 0.0)
edgey = ((n.get_dim('y') == 1.0)|(n.get_dim('y') == 0.0))
labels = np.zeros(n.lennodes)
labels[edgex] = 1; labels[edgey] = 2
labnam = {"mass":0,"edgex":1,"edgey":2}
n.add_labels(labels,labnam)

h = 3

A = make_stiffness(n,2,h,'x',labelname='mass') + make_stiffness(n,2,h,'y',labelname='mass')

A += make_stiffness(n,0,h,'x',labelname='edgex')
A += make_stiffness(n,0,h,'x',labelname='edgey')

print('setup done')

b = np.zeros(n.lennodes)
b[edgex] = 100
b[edgey] = 110

#T = svdsolve(A,b)
from scipy.sparse import linalg
As = sparse.csr_matrix(A)
#T,_ = linalg.bicg(A,b)
#Ts = T.reshape((N,N))
#plt.imshow(Ts,interpolation='bicubic')
#plt.colorbar()
#plt.show()

def func(v):
    prod = As@v
    return prod - b

from scipy.optimize import root

guess = np.zeros(b.shape)

sol = root(func, guess, method='krylov', options={'disp': True},tol=0.000005)

Ts = sol.x.reshape((N,N))
plt.imshow(Ts,interpolation='bicubic')
plt.show()

from pde import Diff, OpList

ops = []

ops.append(Diff(n,'mass','x',2,0.0))
ops.append(Diff(n,'mass','y',2,0.0))
ops.append(Diff(n,'edgex','x',0,100))
ops.append(Diff(n,'edgey','y',0,110))

opl = OpList(ops)
