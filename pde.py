from Node2 import *
from scipy.sparse import csr_matrix
import numpy as np


class Diff:
    islinop = True
    def __init__(self,nodes,label,dim,D,target,h=3):
        self.nodes = nodes
        self.label = label
        self.dim = dim
        self.h = h
        self.D = D
        #self.target = target
        self.mat = make_stiffness(nodes,D,h,dim,labelname=label)
        self.target = np.zeros(nodes.lennodes)
        self.target[nodes.get_label_ids(label)] = target

    def forward(self,vec):
        return self.mat@vec

class LinOp:
    def __init__(self,mat,targ):
        self.mat = mat
        self.targ = targ
    
    def solve(self):
        return np.linalg.solve(self.mat,self.targ)

    def __add__(self,other):
        return LinOp(self.mat + other.mat,self.targ + other.targ)



class OpList:
    def __init__(self,oplist=[]):
        self.oplist = oplist

    def forward(self,vec):
        prod = np.zeros(vec.shape)
        targ = np.zeros(vec.shape)
        for op in self.oplist:
            targ+=op.target
            prod += op.forward(vec)
        return prod, targ

    def combine_linops(self):
        lop = None
        nlop = []
        for op in self.oplist:
            if op.islinop and not lop:
                lop = LinOp(op.mat,op.target)
            if op.islinop and lop:
                lop += LinOp(op.mat,op.target)
            else:
                nlop.append(op)
        self.lop = lop
        self.nlop = nlop

    def solve(self):
        if self.nlop == [] and self.lop:
            return self.lop.solve()



