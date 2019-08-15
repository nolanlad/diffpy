from nodes import Node

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

class Space2:
    def __init__(self,dims,dep_dims):
        self.dimnames = dims
        self.dependent_dims = dep_dims
        self.count = 0

    def nodegen(self,val):
        node = Node(self,self.count,val)
        self.count+=1
        return node
