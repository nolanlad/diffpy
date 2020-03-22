import numpy as np 
import matplotlib.pyplot as plt 
from utilities import factorial
from utilities import *

def argmins(x,N):
    sor = np.argsort(x)
    return sor[:N]

x = np.linspace(0,1,10)
y = np.linspace(0,1,10)

xx,yy = np.meshgrid(x,y)
xx = xx.reshape(-1)
yy = yy.reshape(-1)

x1 = xx[25]
y1 = yy[25]

#xx = xx - x1
#yy = yy - y1

from Node import Nodes

n = Nodes({"x":xx,"y":yy})

def diff(j,nods,D,h):
    N = n_combos(3,2)
    n = nods
    clos = n.get_neighbors2(j,['x','y'],N)

    nod = n.get_node(j)
    '''
    n.plot('x','y','ro')
    clos.plot('x','y','bo')
    plt.show()
    '''
    diffsx = clos.get_dim('x') - nod.get_dim('x')
    diffsy = clos.get_dim('y') - nod.get_dim('y')

    A = np.zeros((N,N))

    for i in range(N):
        A[:,i] = make_combos([diffsx[i],diffsy[i]],2)

    #print(A)
    b = [0,0,factorial(D),0,0,0]
    #return A
    #from diff import svdsolve
    soln = svdsolve(A,b)
    #print(soln)
    return clos.i,soln


i = 25
M = n.lennodes
A = np.zeros((M,M))
for i in range(1,n.lennodes):
    #print(i)
    ids, soln = diff(i,n,1,1)
    A[i][ids] = soln


z = (n.get_dim('x')**2)*n.get_dim('y')
dza = np.matmul(A,z)
ddza = np.matmul(A,dza)
dz = 2*(n.get_dim('x'))*n.get_dim('y')
ddz =  2*n.get_dim('y')

d = np.dot(z[ids],soln)

plt.plot(ddza,ddz,'bo')
plt.show()
'''
clos = n.get_neighbors2(0,['x','y'],6)
n.plot('x','y','bo')
clos.plot('x','y','ro')
plt.show()
'''
#print(d)

#print(dz[i])