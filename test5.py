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

# plt.imshow(Ay)
# plt.show()

z = (n.get_dim('x')**2)*n.get_dim('y')
dzdy = (n.get_dim('x')**2)
dza = np.matmul(A,z)
ddza = np.matmul(A,dza)
dz = 2*(n.get_dim('x')*n.get_dim('y'))
ddz =  2*n.get_dim('y')

d = np.dot(z[ids],soln)

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

