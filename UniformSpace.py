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

        

    
