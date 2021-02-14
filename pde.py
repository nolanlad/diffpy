from Node2 import *
from scipy.sparse import csr_matrix


class Diff:
    islinop = True
    def __init__(self,nodes,label,dim,D,target,h=3):
        self.nodes = nodes
        self.label = label
        self.dim = dim
        self.h = h
        self.D = D
        self.target = target
        self.mat = make_stiffness(nodes,D,h,dim,labelname=label)
        self.target = np.zeros(nodes.lennodes)
        self.target[nodes.get_label_ids(label)] = target

    def forward(self,vec):
        return self.mat@vec




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


