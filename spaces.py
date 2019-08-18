from nodes import Node
from RegularNode import RegularNode

# class Space:
#     def __init__(self):
#         self.dimnames = []
#         self.dependent_dims = []
#         self.count = 0
        
#     def add_dimension(self,name):
#         self.dimnames.append(name)

#     def add_non_fixed_dimension(self,name):
#         self.dependent_dims.addend(name)

#     def nodegen(self,val):
#         node = Node(self,self.count,val)
#         self.count+=1
#         return node

class Space:
    def __init__(self,dims,dep_dims):
        self.dimnames = dims
        self.dependent_dims = dep_dims
        self.count = 0

    def nodegen(self,val):
        node = Node(self,self.count,val)
        self.count+=1
        return node

class RegularSpace(Space):
    def __init__(self,dims,dep_dims):
        Space.__init__(self,dims,dep_dims)

    def nodegen(self,val):
        node = RegularNode(self,self.count,val)
        self.count+=1
        return node

class IrregularSpace(Space):
    def __init__(self,dims,dep_dims):
        Space.__init__(self)

    def nodegen(self,val):
        node = Node(self,self.count,val)
        self.count+=1
        return node