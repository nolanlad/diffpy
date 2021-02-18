import numpy as np 
import matplotlib.pyplot as plt 
from utilities import *
from Node2 import *
from pde import Diff,Diff2,OpList

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

n = Nodes({"x":xx,"y":yy})

edgex = (n.get_dim('x') == 1.0)|(n.get_dim('x') == 0.0)
edgey = ((n.get_dim('y') == 1.0)|(n.get_dim('y') == 0.0))
labels = np.zeros(n.lennodes)
labels[edgex] = 1; labels[edgey] = 2
labnam = {"mass":0,"edgex":1,"edgey":2}
n.add_labels(labels,labnam)

h = 3

ops = []

print("mass 1")
ops.append(Diff2(n,'mass','x',2,0.0))
print('mass 2')
ops.append(Diff2(n,'mass','y',2,0.0))
ops.append(Diff2(n,'edgex','x',0,100))
ops.append(Diff2(n,'edgey','y',0,110))

opl = OpList(ops)

opl.combine_linops()

T = opl.solve()

Ts = T.reshape((N,N))

plt.imshow(Ts,interpolation='bicubic')
plt.show()