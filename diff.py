#from nodes import Node, NodeWalker
from IrregularNode import IrregularNode, IrregularNodeWalker
from RegularNode import RegularNode, RegularNodeWalker
from utilities import *
from spaces import Space, RegularSpace
import numpy as np

def differentiate(node,dim2,N,h):
    ndims = len(node.refspace.dimnames)+1
    dimms = tuple(node.refspace.dimnames)
    p = int(n_combos(ndims,N+h))
    nodes = node.get_neighbors(p,dim2)
    nodes.append(node)
    stiffness = np.zeros((int(p),int(p)))
    for i in range(p):
        dif = nodes[i].get_vals(dimms) - node.get_vals(dimms)
        stiffness[:,i] = make_combos(dif,N+h)
    thelist = make_combos2(node.refspace.dimnames,N+h)
    b = np.zeros(int(p))
    b[thelist.index(dim2*N)] = factorial(N)
    ids = np.array([n.N for n in nodes])
    soln = np.linalg.solve(stiffness,b)
    return ids,soln

def differentiate_debug(node,dim2,N,h):
    ndims = len(node.refspace.dimnames)+1
    dimms = tuple(node.refspace.dimnames)
    p = int(n_combos(ndims,N+h))
    nodes = node.get_neighbors(p,dim2)
    nodes.append(node)
    stiffness = np.zeros((int(p),int(p)))
    diffs = []
    for i in range(p):
        dif = nodes[i].get_vals(dimms) - node.get_vals(dimms)
        diffs.append(dif)
        #print(dif)
        stiffness[:,i] = make_combos(dif,N+h)
    thelist = make_combos2(node.refspace.dimnames,N+h)
    b = np.zeros(int(p))
    b[thelist.index(dim2*N)] = factorial(N)
    ids = np.array([n.N for n in nodes])
    print(thelist)
    #print(stiffness)
    soln = solve_linearly_dependent(stiffness,b)
    return ids,soln
    #return diffs,stiffness

def differentiate_debug2(node,dim2,N,h):
    ndims = len(node.refspace.dimnames)+1
    dimms = tuple(node.refspace.dimnames)
    p = int(n_combos(ndims,N+h))
    nodes = node.get_neighbors(p,dim2)
    nodes.append(node)
    stiffness = np.zeros((int(p),int(p)))
    diffs = []
    for i in range(p):
        dif = nodes[i].get_vals(dimms) - node.get_vals(dimms)
        diffs.append(dif)
        #print(dif)
        stiffness[:,i] = make_combos(dif,N+h)
    thelist = make_combos2(node.refspace.dimnames,N+h)
    b = np.zeros(int(p))
    b[thelist.index(dim2*N)] = factorial(N)
    ids = np.array([n.N for n in nodes])
    print(thelist)
    #print(stiffness)
    #soln = solve_linearly_dependent(stiffness,b)
    #return ids,soln
    return diffs,stiffness


def nodesquare():
    space = Space()
    space.add_dimension("x")
    space.add_dimension("y")
    nodelist = []
    for i in range(-1,2):
        for j in range(-1,2):
            nodelist.append(space.nodegen([i,j]))
    for i in [0,1,2,3,4,5,6,7,8]:
        nodelist[4].link(nodelist[i])
    return nodelist

def nodesquare2():
    space = RegularSpace(('x','y'),('z'))
    #space.add_dimension("x")
    #space.add_dimension("y")
    nodelist = []
    for i in range(-1,2):
        for j in range(-1,2):
            nodelist.append(space.nodegen([i,j]))
    for i in [0,1,2,3,4,5,6,7,8]:
        nodelist[4].link(nodelist[i])
    return nodelist

def node_rect():
    space = RegularSpace(('x','y'),('T'))
    x = np.linspace(0,9,10)
    y = x
    xx,yy = np.meshgrid(x,y)
    nodes = []
    for i in range(10):
        row = []
        for j in range(10):
            node = space.nodegen([xx[i][j],yy[i][j]])
            row.append(node)
        nodes.append(row)
    for i in range(1,9):
        for j in range(1,9):
            nodes[i][j].link(nodes[i][j+1])
            nodes[i][j].link(nodes[i][j-1])
    for i in range(1,9):
        for j in range(1,9):
            nodes[i][j].link(nodes[i-1][j])
            nodes[i][j].link(nodes[i+1][j])
            nodes[i][j].link(nodes[i+1][j+1])
            nodes[i][j].link(nodes[i-1][j-1])
            nodes[i][j].link(nodes[i-1][j+1])
            nodes[i][j].link(nodes[i+1][j-1])
    #nodes2 = []
    #for r in nodes:
    #    nodes2.extend(r)
    return nodes

def nodelinspace():
    nodes = []
    space = Space2(('x','y'),('T'))
    for i in np.arange(10):
        node = space.nodegen([i,np.random.random()])
        nodes.append(node)
    for i in range(1,9):
        nodes[i].link(nodes[i-1])
        nodes[i].link(nodes[i+1])
    return nodes

def nodelinspace2():
    nodes = []
    space = Space2(('x'),('T'))
    for i in np.arange(20):
        node = space.nodegen([i])
        nodes.append(node)
    for i in range(1,19):
        nodes[i].link(nodes[i-1])
        nodes[i].link(nodes[i+1])
    return nodes
    
    
def nodelinspace3():
    nodes = []
    space = RegularSpace(('x','y'),('T'))
    for i in np.arange(20):
        node = space.nodegen([i,0])
        nodes.append(node)
    for i in range(1,19):
        nodes[i].link(nodes[i-1])
        nodes[i].link(nodes[i+1])
    return nodes