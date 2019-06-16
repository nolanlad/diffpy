from UniformGraph import UniformGraph1D
import numpy as np
import matplotlib.pyplot as plt

def first_derivative_matrix(graph):
    mat = np.zeros((graph.N,graph.N))
    for i in range(graph.N):
        mat[i] = graph.first_derivative(i)
    return mat

def second_derivative_matrix(graph):
    mat = np.zeros((graph.N,graph.N))
    for i in range(graph.N):
        mat[i] = graph.second_derivative(i)
    return mat

def zero_derivative_matrix(graph):
    mat = np.eye(graph.N)
    return mat

def example1():
    '''
    Solve the equation

    x' + x = 0, x'(1) = 1
    '''
    space = UniformGraph1D(0,1,1000)
    b = np.zeros(1000)
    # x' + x = 0
    stiff = first_derivative_matrix(space) + zero_derivative_matrix(space)
    #x'(1) = 1
    stiff[-1] = space.first_derivative(999)
    b[-1] = 1
    #solve
    c = np.linalg.solve(stiff,b)
    return c

def example2():
    '''
    Solve the non transient heat equation
    with heat generation

    a d2T/dx2 = Qgen, T(0) = T(1) = 0 C
    '''
    space = UniformGraph1D(0,1.0,100)
    b = np.ones(100)*-0.1 #qgen 
    alpha = 1.27e-4
    stiff = alpha*second_derivative_matrix(space)
    stiff[0] = 0
    stiff[-1] = 0
    stiff[0][0] = 1
    stiff[-1][-1] = 1
    b[0] = 273
    b[-1] = 273
    c = np.linalg.solve(stiff,b)
    return c

def example3():
    '''
    solve the ode:
    x'' + x = 0, x(0) = 1, x'(2*pi) = 1

    anylitical solution:
    x(t) = sin(t) + 1
    '''
    space = UniformGraph1D(0,2*np.pi,100)
    A = second_derivative_matrix(space) + zero_derivative_matrix(space)
    b = np.ones(100)
    A[0] = 0
    A[-1] = 0
    A[0][0] = 1
    A[-1] = space.first_derivative(99)
    soln = np.linalg.solve(A,b)
    plt.plot(soln,'bo')
    plt.show()
    
