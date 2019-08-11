from nodes import Node, NodeWalker
from utilities import *
from spaces import Space
import numpy as np

def differentiate(node,dim2,N,h):
    ndims = len(node.refspace.dimnames)+1
    dimms = tuple(node.refspace.dimnames)
    p = int(n_combos(ndims,N+h))
    walker = NodeWalker(node,p)
    walker.add_variance_var(dim2)
    diffs = []
    nodes = [n for n in walker]
    nodes.append(node)
    for n in nodes:
        diffs.append(n.get_vals(dimms) - node.get_vals(dimms) )
    stiffness = np.zeros((int(p),int(p)))
    for i in range(p):
        stiffness[:,i] = make_combos(diffs[i],N+h)
    thelist = make_combos2(node.refspace.dimnames,N+h)
    b = np.zeros(int(p))
    b[thelist.index(dim2*N)] = factorial(N)
    ids = [n.N for n in nodes]
    soln = np.linalg.solve(stiffness,b)
    return ids,soln

def nodesquare():
    space = Space()
    space.add_dimension("x")
    space.add_dimension("y")
    #space.add_dimension("z")
    nodelist = []
    for i in range(-1,2):
        for j in range(-1,2):
            nodelist.append(space.nodegen([i,j]))
    for i in [0,1,2,3,4,5,6,7,8]:
        nodelist[4].link(nodelist[i])
    return nodelist
