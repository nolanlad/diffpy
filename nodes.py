import numpy as np

class Node:
    def __init__(self,refspace,N,val):
        self.refspace = refspace
        self.N = N
        self.val = np.array(val)
        self.neighbors = []
        
    def get_val(self,dimname):
        ind = self.refspace.dimnames.index(dimname)
        return self.val[ind]

    def get_vals(self,dims):
        return np.array([self.get_val(d) for d in dims])

    def link(self,other):
        self.neighbors.append(other)

    def get_neighbors(self,N,var_dim):
        walker = NodeWalker(self,N)
        walker.add_invariance(var_dim)
        nodes = [n for n in walker]
        return nodes

class NodeWalker:
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
        self.count =0

    def __iter__(self):
        return self

    def __next__(self):
        '''
        Please someone who knows graph theory make this not suck
        '''
        if self.count == self.N-1:
            raise StopIteration()
        for n in self.layer:
            for n2 in n.neighbors:
                if n2.N not in self.visited:
                    if self.variance_dim != None and (n2.get_val(self.variance_dim) -
                        self.startnode.get_val(self.variance_dim)) != 0:
                        self.visited.append(n2.N)
                        self.count+=1
                        return n2
                    else:
                        self.visited.append(n2.N)
        newlayer = []
        for l in self.layer:
            newlayer.extend(l.neighbors)
        self.layer = newlayer
        
        for n in self.layer:
            for n2 in n.neighbors:
                if n2.N not in self.visited:
                    if (n2.get_val(self.variance_dim) - self.startnode.get_val(self.variance_dim)) != 0:
                        self.visited.append(n2.N)
                        self.count+=1
                        return n2
                    else:
                        self.visited.append(n2.N)
                        
    def next(self):
        return self.__next__()

    def add_variance_var(self,dim):
        self.variance_dim = dim
