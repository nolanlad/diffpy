import numpy as np
import matplotlib.pyplot as plt

class UniformGraph1D:
    def __init__(self,x0,xf,N):
        self.x0 = x0
        self.xf = xf
        self.N = N
        self.gran = (xf-x0)/N
        
    def getx(self,ind):
        return ind*gran

    def needs_forward(self,ind,degree=2):
        return ind == 0

    def needs_backward(self,ind,degree=2):
        return ind == (self.N-1)
    
    def central_derivative(self,ind,degree=2):
        out = np.zeros(self.N)
        out[ind - 1] = (-1/2)*(1/self.gran)
        out[ind + 1] = (1/2)*(1/self.gran)
        return out

    def forward_derivative(self,ind,degree=2):
        out = np.zeros(self.N)
        out[ind] = (-3/2)*(1/self.gran)
        out[ind+1] = 2*(1/self.gran)
        out[ind+2] = (-1/2)*(1/self.gran)
        return out

    def backward_derivative(self,ind,degree=2):
        out = np.zeros(self.N)
        out[ind] = (3/2)*(1/self.gran)
        out[ind-1] = -2*(1/self.gran)
        out[ind-2] = (1/2)*(1/self.gran)
        return out

    def central_derivative2(self,ind,degree=2):
        out = np.zeros(self.N)
        out[ind-1:ind+1+1] = np.array([
            1,-2,1])/(self.gran**2)
        return out

    def backward_derivative2(self,ind,degree=2):
        out = np.zeros(self.N)
        out[ind-3:ind+1] = np.array([
            -1,4,-5,2])/(self.gran**2)
        return out

    def forward_derivative2(self,ind,degree=2):
        out = np.zeros(self.N)
        out[ind:ind+4] = np.array([
            2,-5,4,-1])/(self.gran**2)
        return out
    
    def first_derivative(self,ind,degree=2):
        if self.needs_forward(ind):
            return self.forward_derivative(ind,degree=degree)
        if self.needs_backward(ind):
            return self.backward_derivative(ind,degree=degree)
        return self.central_derivative(ind,degree=degree)

    def second_derivative(self,ind,degree=2):
        if self.needs_forward(ind):
            return self.forward_derivative2(ind,degree=degree)
        if self.needs_backward(ind):
            return self.backward_derivative2(ind,degree=degree)
        return self.central_derivative2(ind,degree=degree)

    




