import numpy as np
import matplotlib.pyplot as plt

class Node:
    def __init__(self,refspace,N,val):
        self.refspace = refspace
        self.N = N
        self.val = np.array(val)
        self.neighbors = []
        self.is_edge = False
        
    def get_val(self,dimname):
        ind = self.refspace.dimnames.index(dimname)
        return self.val[ind]

    def get_vals(self,dims):
        return np.array([self.get_val(d) for d in dims])

    def link(self,other):
        if other not in self.neighbors:
            self.neighbors.append(other)
        if self not in other.neighbors:
            other.neighbors.append(self)


class NodeWalker:
    '''
    tbh I hate this solution to this problem
    it's not elegant or efficient. Please someone
    who is actually good at computer science make
    a class that isn't a steaming pile of shit
    '''
    def __init__(self,node):
        self.startnode = node
        self.visited = [node.N]
        self.layer = [node]
        #self.N = N
        self.variance_dim = None
        self.count =1

    def __iter__(self):
        return self

    def __next__(self):
        '''
        Please someone who knows graph theory make this not suck
        '''
        for n in self.layer:
            for n2 in n.neighbors:
                if n2.N not in self.visited:
                    self.visited.append(n2.N)
                    self.count+=1
                    return n2
        newlayer = []
        for l in self.layer:
            newlayer.extend(l.neighbors)
        self.layer = newlayer
        
        for n in self.layer:
            for n2 in n.neighbors:
                if n2.N not in self.visited:
                    self.visited.append(n2.N)
                    self.count+=1
                    return n2
                        
    def next(self):
        return self.__next__()

    def add_variance_var(self,dim):
        self.variance_dim = dim



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
    def get_neighbors(self,i,dims,N):
        mag = np.zeros(self.size())
        nod = self.get_node(i)
        diffsx = self.get_dim('x') - nod.get_dim('x')
        diffsy = self.get_dim('y') - nod.get_dim('y')
        
        for d in dims:
            dum = self.get_dim(d)
            dum2 = nod.get_dim(d)
            mag = mag + (dum-dum2)**2
        return  SubNode(np.argsort(np.sqrt(mag))[:N],self)
    def get_neighbors2(self,i,dims,N):
        mag = np.zeros(self.size())
        nod = self.get_node(i)
        diffsx = self.get_dim('x') - nod.get_dim('x')
        diffsy = self.get_dim('y') - nod.get_dim('y')
        criterion = (np.abs(diffsx)+1)*np.arctan(0.1 + np.abs(diffsy))
        # plt.imshow(criterion.reshape((10,10)))
        # plt.colorbar()
        # plt.show()
        return  SubNode(np.argsort(criterion)[:N],self)
        
        # ang = np.argsort(np.arctan(np.abs(diffsy/diffsx)))
        # line = np.argsort(np.abs(diffs[ang]))
        

    def plot(self,dim1,dim2,*args,**kwargs):
        plt.plot(self.get_dim(dim1),self.get_dim(dim2),*args,**kwargs)
'''
    def diff(self,i):
        A = np.zeros(6,6)
        clos = self.get_neighbors(i,['x','y'],6)
        nod = self.get_node(i)
        diffsx = clos.get_dim('x') - nod.get_dim('x')
        diffsx = clos.get_dim('y') - nod.get_dim('y')
        for i in range(6):
            return
'''



class SubNode(Nodes):
    def __init__(self,i,refspace):
        self.refspace = refspace
        self.i = i
    def get_dim(self,name):
        return self.refspace.get_dim(name)[self.i]