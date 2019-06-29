import numpy as np


'''============================================================
Helper functions to find the coefficents of finite differencing
============================================================'''

def factorial(N):
    if N == 0:
        return 1
    tot = 1
    for i in range(1,N+1):
        tot*=i
    return tot

def get_p(N,h):
    p = ((2*((N+1)/2))-2 + h)/2
    return int(p)

def get_central_diff_coeffs(N,h):
    p = ((2*((N+1)/2))-2 + h)/2
    p = int(p)
    A = np.zeros((2*p+1,2*p+1))
    for i in range(-p,p+1):
        for j in range(2*p+1):
            A[j][i+p] = (i**j)
    b = np.zeros(2*p+1)
    b[N] = factorial(N)
    coeffs = np.linalg.solve(A,b)
    return coeffs

def get_forward_diff_coeffs(N,h):
    p = int(N+h)
    A = np.zeros((p,p))
    for i in range(0,p):
        for j in range(0,p):
            A[j][i] = (i**j)
    b = np.zeros(p)
    b[N] = factorial(N)
    coeffs = np.linalg.solve(A,b)
    return coeffs

def get_backward_diff_coeffs(N,h):
    coeffs = get_forward_diff_coeffs(N,h)
    return coeffs[::-1]
'''
================================================================
'''


class UniformSpace1D:
    def __init__(self,x0,xf,N):
        self.x0 = x0
        self.xf = xf
        self.N = N
        self.gran = (xf-x0)/N

    def get_ind(ind):
        return ind

    def needs_forward(self,ind,N=1,h=2):
        return ind-get_p(N,h) < 0

    def needs_backward(self,ind,N=1,h=2):
        return ind+get_p(N,h) > (self.N-1)

    def central_diff(self,A,ind,N,h):
        coeffs = get_central_diff_coeffs(N,h)
        p = len(coeffs)//2
        for i in range(-p, p+1):
            A[ind][ind+i] = coeffs[i+p]
            
    def front_diff(self,A,ind,N,h):
        coeffs = get_forward_diff_coeffs(N,h)
        p = len(coeffs)
        for i in range(p):
            A[ind][ind+i] = coeffs[i]

    def back_diff(self,A,ind,N,h):
        coeffs = get_backward_diff_coeffs(N,h)
        p = len(coeffs)
        for i in range(p):
            A[ind][ind-i] = coeffs[-1-i]

    def diff(self,A,ind,N,h):
        if self.needs_forward(ind,N=N,h=h):
            self.front_diff(A,ind,N,h)
        elif self.needs_backward(ind,N=N,h=h):
            self.back_diff(A,ind,N,h)
        else:
            self.central_diff(A,ind,N,h)
'''
class Node:
    def __init__(self,loc,is_edge=False):
        self.loc = loc
        self.neighbors = []
        self.is_edge = is_edge
    def setloc(self,loc):
        self.loc = np.array(loc)
'''

class Space:
    def __init__(self):
        self.dimnames = []
        self.count = 0
        
    def add_dimension(self,name):
        self.dimnames.append(name)

    def nodegen(self,val):
        node = Node(self,self.count,val)
        self.count+=1
        return node


class Node:
    def __init__(self,refspace,N,val):
        self.refspace = refspace
        self.N = N
        self.val = np.array(val)
        self.neighbors = []
        
    def get_val(self,dimname):
        ind = self.refspace.dimnames.index(dimname)
        return self.val[ind]

    def diff(self,dim1,dim2,N,h):
        p = N+h
        A = np.zeros((p,p))
        b = np.zeros(p)
        A[0][0] = 1 #all f_i cancel to zero
        current_node = self
        indexed_nodes = [self.N]
        walker = NodeWalker(self,p)
        column = 1
        y = [self.get_val(dim1)]
        for n in walker:
            h = n.get_val(dim2) - self.get_val(dim2)
            y.append(n.get_val(dim1))
            for i in range(p):
                A[i][column] = h**i
            column+=1
        b[N] = factorial(N)
        coeffs = np.linalg.solve(A,b)
        return coeffs, np.array(walker.visited)
            
            
            

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

    def __iter__(self):
        return self

    def next(self):
        '''
        Please someone who knows graph theory make this not suck
        '''
        if len(self.visited) == self.N:
            raise StopIteration()
        for n in self.layer:
            for n2 in n.neighbors:
                if n2.N not in self.visited:
                    self.visited.append(n2.N)
                    return n2
        newlayer = []
        for l in self.layer:
            newlayer.extend(l.neighbors)
        self.layer = newlayer
        
        for n in self.layer:
            for n2 in n.neighbors:
                if n2.N not in self.visited:
                    self.visited.append(n2.N)
                    return n2

                
def nodelinspace(start,stop,N):
    '''Delete this'''
    space = Space()
    space.add_dimension("x")
    space.add_dimension("y")
    node_list = [space.nodegen([i,i**2]) for i in np.linspace(start,stop,N)]
    node_list[0].neighbors.append(node_list[1])
    node_list[-1].neighbors.append(node_list[-2])
    for i in range(1,len(node_list)-1):
        node_list[i].neighbors.append(node_list[i-1])
        node_list[i].neighbors.append(node_list[i+1])
    return node_list
        
            
        
        
    
    
