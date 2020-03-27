import numpy as np
import matplotlib.pyplot as plt
from utilities import *


class Nodes:
    def __init__(self,space):
        self.dims = np.array([k for k in space])
        self.lendims = len(self.dims)
        self.lennodes = len(space[self.dims[0]]) 
        self.idims = dict(zip(self.dims,range(self.lendims)))
        self.space = np.zeros((self.lendims,self.lennodes))
        for i,k in enumerate(self.dims):
            self.space[i] = space[k]

    def add_labels(self,labs,labnames):
        '''labs: labels for nodes
        labnames: dictionary containing meanings of labels
        '''
        self.labels = labs
        self.labnames = labnames
    def size(self):
        return self.lennodes
        #return len(self.space[self.dims[0]])
    def get_dim(self,name):
        return self.space[self.idims[name]]
    def add_dims(self,dim):
        'BAD! WIP'
        for k in dim.keys():
            self.space[k] = dim[k]
    def get_node(self,i):
        return SubNode(i,self)

    def get_neighbors_x(self,i,dims,N):
        mag = np.zeros(self.size())
        nod = self.get_node(i)
        diffsx = self.get_dim('x') - nod.get_dim('x')
        diffsy = self.get_dim('y') - nod.get_dim('y')
        criterion = (np.abs(diffsx)+1)*np.arctan(0.1 + np.abs(diffsy))
        return  SubNode(np.argsort(criterion)[:N],self)

    def get_neighbors_y(self,i,N):
        mag = np.zeros(self.size())
        nod = self.get_node(i)
        diffsx = self.get_dim('x') - nod.get_dim('x')
        diffsy = self.get_dim('y') - nod.get_dim('y')
        criterion = (np.abs(diffsy)+1)*np.arctan(0.1 + np.abs(diffsx))
        return  SubNode(np.argsort(criterion)[:N],self)
    
    def get_neighbors(self,dim,i,N):
        nod = self.get_node(i)
        diffvariant = self.get_dim(dim) - nod.get_dim(dim)
        criterion = (np.abs(diffvariant)+1)
        invariantdims = [k for k in self.dims if k!=dim]
        for idim in invariantdims:
            dif = self.get_dim(idim) - nod.get_dim(idim)
            criterion = criterion*np.arctan(0.1 + np.abs(dif))
        return  SubNode(np.argsort(criterion)[:N],self)

    def plot(self,dim1,dim2,*args,**kwargs):
        for l in set(self.labels):
            w = (self.labels == l)
            plt.plot(self.get_dim(dim1)[w],self.get_dim(dim2)[w],'o',*args,**kwargs)
        plt.legend([l for l in self.labnames])

def diff_x(j,nods,D,h):
    N = n_combos(3,2)
    n = nods
    #clos = n.get_neighbors_x(j,['x','y'],N)
    clos = n.get_neighbors('x',j,N)

    nod = n.get_node(j)
    diffsx = clos.get_dim('x') - nod.get_dim('x')
    diffsy = clos.get_dim('y') - nod.get_dim('y')

    A = np.zeros((N,N))

    for i in range(N):
        A[:,i] = make_combos([diffsx[i],diffsy[i]],2)

    thelist = make_combos2('xy',D+h)
    b = np.zeros(int(N))
    b[thelist.index('x'*D)] = factorial(D)
    soln = svdsolve(A,b)
    return clos.i,soln

def diff_gen(j,nods,D,h,dim):
    N = n_combos(3,D+h)
    n = nods
    #clos = n.get_neighbors_x(j,['x','y'],N)
    clos = n.get_neighbors(dim,j,N)

    nod = n.get_node(j)
    diffsx = clos.get_dim('x') - nod.get_dim('x')
    diffsy = clos.get_dim('y') - nod.get_dim('y')

    A = np.zeros((N,N))

    for i in range(N):
        A[:,i] = make_combos([diffsx[i],diffsy[i]],D+h)

    thelist = make_combos2('xy',D+h)
    b = np.zeros(int(N))
    b[thelist.index(dim*D)] = factorial(D)
    soln = svdsolve(A,b)
    return clos.i,soln


def diff_y(j,nods,D,h):
    N = n_combos(3,D+h)
    n = nods
    clos = n.get_neighbors_y(j,N)

    nod = n.get_node(j)
 
    diffsx = clos.get_dim('x') - nod.get_dim('x')
    diffsy = clos.get_dim('y') - nod.get_dim('y')

    A = np.zeros((N,N))

    for i in range(N):
        A[:,i] = make_combos([diffsx[i],diffsy[i]],D+h)
    thelist = make_combos2('xy',D+h)
    b = np.zeros(int(N))
    b[thelist.index('y'*D)] = factorial(D)

    soln = svdsolve(A,b)
    return clos.i,soln

def make_stiffness(nodes,D,h,dim,labelname = None):
    n = nodes
    M = n.lennodes
    if D == 0:
        return np.eye(M)
    A = np.zeros((M,M))
    if not labelname:
        for i in range(0,n.lennodes):
            ids, soln = diff_gen(i,n,D,h,dim)
            A[i][ids] = soln
    else:
        ind = n.labnames[labelname]
        w = np.where(n.labels == ind)[0]
        for i in w:
            ids, soln = diff_gen(i,n,D,h,dim)
            A[i][ids] = soln

    return A

class SubNode(Nodes):
    def __init__(self,i,refspace):
        self.refspace = refspace
        self.i = i
    def get_dim(self,name):
        return self.refspace.get_dim(name)[self.i]