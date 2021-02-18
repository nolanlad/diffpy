import numpy as np 
import matplotlib.pyplot as plt 
from utilities import *
from Node2 import *
from Node3 import *
from pde import Diff,OpList

'''
Solve for steady state solution for heat conduction in a
square of material with adiabatic walls and constant temperatures
at top and bottom
'''
N=31
x = np.linspace(0,2,N)
y = np.linspace(0,2,N)

xx,yy = np.meshgrid(x,y)
xx = xx.reshape(-1)
yy = yy.reshape(-1)

nods = Nodes({"x":xx,"y":yy})



sol,inds,A,b = test(nods.get_node(41),nods,'x',2,0)
sol,inds,A,b = test(nods.get_node(17),nods,'x',2,0)

A = make_stiffness_general(nods,1,3,'x')

res = -1*A@nods.get_dim('x')

plt.imshow(res.reshape((N,N)))
plt.colorbar()
plt.show()

zz = xx*xx*yy

dza = -1*A@zz
dz = 2*xx*yy

plt.imshow(dza.reshape((N,N)))
plt.colorbar()
plt.show()

A = make_stiffness_general(nods,1,3,'y')
#A2 = make_stiffness(nods,,1,'y')

res2 = -1*A@nods.get_dim('y')

plt.imshow(res2.reshape((N,N)))
plt.colorbar()
plt.show()

zz = xx*xx*yy

dza = -1*A@zz
dz = xx*xx

plt.imshow(dza.reshape((N,N)))
plt.colorbar()
plt.show()



'''
n = nods.get_node(N)
soln,inds = get_node_soln(n,nods,'x',1,0)

x = nods.get_dim('x')
y = nods.get_dim('y')
z = x*x*y*y
y2 = nods.get_dim('y')[N]
x2 = nods.get_dim('x')[N]

Z = make_stiffness_5(nods,2,3,'x')

ddzdx_a = Z@z

ddzdx = 2*x*y*y

plt.imshow(ddzdx.reshape((25,25))-ddzdx_a.reshape((25,25)))
plt.colorbar()
plt.show()
plt.imshow(ddzdx_a.reshape((25,25)))
plt.show()
'''
