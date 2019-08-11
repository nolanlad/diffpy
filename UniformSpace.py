import numpy as np
from scipy import linalg
import itertools


'''============================================================
Helper functions to find the coefficents of finite differencing
============================================================'''

def solve_non_square(A,b):
    L,U = linalg.lu(A,permute_l = True)
    trivial_col = 0
    for i in range(len(U)):
        if U[i][i] == 0:
            trivial_col = i
            break
    Ut = np.transpose(U)
    useful_rows = []
    for i in range(U.shape[0]):
        for j in range(U.shape[1]):
            if U[i][j] != 0:
                useful_rows.append(j)
                break
    useful_rows = np.array(useful_rows)
    Ut2 = Ut[useful_rows]
    U2 = np.transpose(Ut2)
    b2 = np.linalg.inv(L)@b
    solnt = np.linalg.solve( U2,b2 )
    soln = np.zeros(U.shape[1])
    soln[useful_rows] = solnt
    return soln
    

def factorial(N):
    if N == 0:
        return 1
    tot = 1
    for i in range(1,N+1):
        tot*=i
    return tot

def n_combos(n,r):
    return factorial(n+r-1)/((factorial(r)*factorial(n-1)))

def make_combos(shit,n):
    shit = list(shit)
    shit.extend([1])
    return [multiply_els(c) for c in \
            itertools.combinations_with_replacement(shit,n)]

def make_combos2(shit2,n):
    shit = list(shit2)
    shit.append(' ')
    return [sum_els(c).strip() for c in \
            itertools.combinations_with_replacement(shit,n)]
    

def multiply_els(iterr):
    acc = 1
    for i in iterr:
        acc*=i
    return acc

def sum_els(iterr):
    acc = ''
    for i in iterr:
        acc+=i
    return acc


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
        self.dependent_dims = []
        self.count = 0
        
    def add_dimension(self,name):
        self.dimnames.append(name)

    def add_non_fixed_dimension(self,name):
        self.dependent_dims.addend(name)

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

    def get_vals(self,dims):
        return np.array([self.get_val(d) for d in dims])

    def link(self,other):
        self.neighbors.append(other)

    def get_neighbors(self,N,var_dim):
        walker = NodeWalker(self,N)
        walker.add_invariance(var_dim)
        nodes = [n for n in walker]
        return nodes
    
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
        self.variance_dim = None
        self.count =0

    def __iter__(self):
        return self

    def __next__(self):
        '''
        Please someone who knows graph theory make this not suck
        '''
        #if len(self.visited) == self.N:
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
                        print(n2.N)
                        #return n2
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
                        print('b')
                        #return n2
                        
    def next(self):
        return self.__next__()

    def add_variance_var(self,dim):
        self.variance_dim = dim

        
        

                
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

def nodesquare2():
    space = Space()
    space.add_dimension("x")
    space.add_dimension("y")
    #space.add_dimension("z")
    nodelist = []
    for i in range(-1,2):
        for j in range(-1,2):
            nodelist.append(space.nodegen([i,j]))
    for i in [1,3,5,7]:
        nodelist[4].link(nodelist[i])
    return nodelist       
    

def test1(node,dim1,dim2,N,h):
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
        stiffness[:,i] = make_combos(diffs[i],2)
    #stiffness[p-1][p-1] = 1
    print(make_combos2(node.refspace.dimnames,N+h))
    print([n.N for n in nodes])
    return stiffness

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

def test1(node,dim1,dim2,N,h):
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
        stiffness[:,i] = make_combos(diffs[i],2)
    #stiffness[p-1][p-1] = 1
    print(make_combos2(node.refspace.dimnames,N+h))
    print([n.N for n in nodes])
    return stiffness


'''
def test2(node,dim1,dim2,N,h):
    ndims = len(node.refspace.dimnames)+1
    dimms = tuple(node.refspace.dimnames)
    p = int(n_combos(ndims,N+h))
    print(p)
    walker = NodeWalker(node,p)
    #walker.add_variance_var(dim2)
    diffs = []
    nodes = [n for n in walker]
    for n in nodes:
        diffs.append(n.get_vals(dimms) - node.get_vals(dimms) )
    stiffness = np.zeros((int(p),int(p)))
    for i in range(p-1):
        stiffness[:,i] = make_combos(diffs[i],1)
    stiffness[0][p-1] = 1
    print(make_combos2(node.refspace.dimnames,N+h))
    print([n.id for n in nodes])
    return stiffness
'''
