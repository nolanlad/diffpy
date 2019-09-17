from Node import NodeWalker

class NodeWalker2:
    def __init__(self,node,var):
        self.origin = node
        self.var = var
        self.walker = NodeWalker(node)
        self.max = 0.0
        self.min = 0.0
    
    def __iter__(self):
        return self
    
    def __next__(self):
        while True:
            node = self.walker.next()

            if not node:
                raise StopIteration()

            diff = self.origin.get_val(self.var) - \
            node.get_val(self.var)

            if diff > self.max:
                self.max = diff
                return node
            if diff < self.min:
                self.min = diff
                return node
        
