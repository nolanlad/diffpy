import numpy as np
from scipy import linalg
import itertools

def solve_non_square(A,b):
    '''Gets one solution to a matrix where m > n'''
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
    b2 = np.matmul(np.linalg.inv(L),b)
    solnt = np.linalg.solve( U2,b2 )
    soln = np.zeros(U.shape[1])
    soln[useful_rows] = solnt
    return soln

def solve_linearly_dependent(matrix,b):
    dept_rows = []
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[0]):
            if i != j:
                inner_product = np.inner(matrix[i],matrix[j])
                norm_i = np.linalg.norm(matrix[i])
                norm_j = np.linalg.norm(matrix[j])
                if np.abs(inner_product - norm_j * norm_i) < 1E-5:
                    dept_rows.extend([i,j])
    idxs = list(set(dept_rows))
    idx = idxs[0]
    matrix2 = np.concatenate((matrix[:idx],matrix[idx+1:]))
    b2 = np.concatenate((b[:idx],b[idx+1:]))
    soln = solve_non_square(matrix2,b2)
    print(np.allclose(np.matmul(matrix,soln),b))
    return soln

def get_linearly_dependent_rows(matrix,b):
    dept_rows = []
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[0]):
            if i != j:
                inner_product = np.inner(matrix[i],matrix[j])
                norm_i = np.linalg.norm(matrix[i])
                norm_j = np.linalg.norm(matrix[j])
                if np.abs(inner_product - norm_j * norm_i) < 1E-5:
                    dept_rows.extend([i,j])
    idxs = list(set(dept_rows))
    return idxs
    

def svdsolve(A,b):
    Ainv = np.linalg.pinv(A)
    return np.matmul(Ainv,b)

def factorial(N):
    if N == 0:
        return 1
    tot = 1
    for i in range(1,N+1):
        tot*=i
    return tot

def n_combos(n,r):
    '''returns number of combinations with replacement of a set '''
    return int(factorial(n+r-1)/((factorial(r)*factorial(n-1))))


def make_combos(shit,n):
    ''' returns combinations of set of floats
    and multiplies them together'''
    shit = list(shit)
    shit.extend([1])
    return [multiply_els(c) for c in \
            itertools.combinations_with_replacement(shit,n)]


def make_combos2(shit2,n):
    '''makes combinations of variables UP TO n'''
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
