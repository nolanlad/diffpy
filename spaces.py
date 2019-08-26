from IrregularNode import IrregularNode
from RegularNode import RegularNode

class Space:
    def __init__(self,dims,dep_dims):
        self.dimnames = dims
        self.dependent_dims = dep_dims
        self.count = 0


class RegularSpace(Space):
    def __init__(self,dims,dep_dims):
        Space.__init__(self,dims,dep_dims)

    def nodegen(self,val):
        node = RegularNode(self,self.count,val)
        self.count+=1
        return node

class IrregularSpace(Space):
    def __init__(self,dims,dep_dims):
        Space.__init__(self,dims,dep_dims)

    def nodegen(self,val):
        node = IrregularNode(self,self.count,val)
        self.count+=1
        return node