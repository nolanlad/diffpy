import numpy as np
from utilities import *
from Node import Node
from matrix import solve_matrix_inf

class RegularNode(Node):
    def __init__(self,refspace,N,val):
        Node.__init__(self,refspace,N,val)

    def get_neighbors(self,N,var_dim):
        walker = RegularNodeWalker(self,N)
        walker.add_variance_var(var_dim)
        nodes = [n for n in walker]
        return nodes
    
    def diff(self,dim,N,h):
        ndims = 1
        p = int(n_combos(ndims+1,N+h))
        nodds = self.get_neighbors(p,dim)
        nodds.append(self)
        stiffness = np.zeros((int(p),int(p)))
        diffs = []
        for i in range(p):
            dif = nodds[i].get_vals([dim]) - self.get_vals([dim])
            diffs.append(dif)
            #print(dif)
            stiffness[:,i] = make_combos(dif,N+h)
        thelist = make_combos2([dim],N+h)
        b = np.zeros(int(p))
        b[thelist.index(dim*N)] = factorial(N)
        ids = np.array([n.N for n in nodds])
        #print(thelist)
        #soln = np.linalg.solve(stiffness,b)
        soln = solve_matrix_inf(stiffness,b)
        return ids,soln

class RegularNodeWalker:
    '''
    tbh I hate this solution to this problem
    it's not elegant or efficient. Please someone
    who is actually good at computer science make
    a class that isn't a steaming pile of shit
    '''
    def __init__(self,node,N):
        self.startnode = node
        self.visited = [node.N]
        self.layer = [node]
        self.N = N
        self.variance_dim = None
        self.count =1

    def add_variance_var(self,dim):
        self.variance_dim = dim
        all_dims = np.array(self.startnode.refspace.dimnames)
        if len(all_dims) == 1:
            self.invariance_dims = None
        else:
            self.invariance_dims = all_dims[np.where(all_dims != dim)]


    def __iter__(self):
        return self
        
    def __next__(self):
        '''
        Please someone who knows graph theory make this not suck
        '''
        if self.count == self.N:
            raise StopIteration()
        for n in self.layer:
            for n2 in n.neighbors:
                if n2.N not in self.visited:
                    if self.invariance_dims != None:
                        if np.all((n2.get_vals(self.invariance_dims) - \
                            self.startnode.get_vals(self.invariance_dims)) == 0):
                            self.visited.append(n2.N)
                            self.count+=1
                            return n2
                        else:
                            self.visited.append(n2.N)
                    else:
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
                    if self.invariance_dims != None:
                        if np.all((n2.get_vals(self.invariance_dims) - \
                            self.startnode.get_vals(self.invariance_dims)) == 0):
                            self.visited.append(n2.N)
                            self.count+=1
                            return n2
                        else:
                            self.visited.append(n2.N)
                    else:
                        self.visited.append(n2.N)
                        self.count+=1
                        return n2