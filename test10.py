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


b = np.zeros(n.lennodes)
b[edgex] = 100
b[edgey] = 400






T = svdsolve(A,b)
Ts = T.reshape((N,N))
plt.imshow(Ts,interpolation='bicubic')
plt.colorbar()
plt.show()
