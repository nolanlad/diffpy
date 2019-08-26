import numpy as np
from scipy import linalg

def is_invertible(A):
    return np.linalg.det(A) == 0

def how_many_solutions(A,b):
    L,U = linalg.lu(A,permute_l = True)
    b2 = np.matmul(np.linalg.inv(L),b)
    (m,n) = A.shape
    if m==n:
        for i in range(m):
            if np.all(U[i] == 0):
                if b2[i] == 0:
                    print('infinite solutions')
                if b2[i] != 0:
                    print('zero solutions')
                
            
        

def solve_infinite_solutions(A,b):
    L,U = linalg.lu(A,permute_l = True)
    b2 = np.matmul(np.linalg.inv(L),b)
    (m,n) = A.shape
    useless_rows = []
    for i in range(m):
        if np.all(U[i] == 0) and b2[i] == 0:
            useless_rows.append(i)
    offset = 0
    for ind in useless_rows:
        U = remove_row_col(U,ind-offset)
        b2 = np.concatenate((b2[:ind -offset],b2[ind+1-offset:]))
        offset+=1
    
    return U,b2

def remove_row_col(A,i):
    A2 = np.concatenate((A[:i],A[i+1:]))
    A2 = np.concatenate((A2[:,:i],A2[:,(i+1):]),axis=1)
    return A2

def remove_row(A,i):
    A2 = np.concatenate((A[:i],A[i+1:]))
    #A2 = np.concatenate((A2[:,:i],A2[:,(i+1):]),axis=1)
    return A2

def remove_col(A,i):
    #A2 = np.concatenate((A[:i],A[i+1:]))
    A2 = np.concatenate((A[:,:i],A[:,(i+1):]),axis=1)
    return A2

def linear_dept_cols(A):
    L,U = linalg.lu(A,permute_l = True)
    Ub = np.zeros(U.shape)
    for i in range(len(U)):
        for j in range(len(U)):
            if np.abs(U[i][j]) >= 1e-9:
                Ub[i][j] = 1
                break
    cols = np.sum(Ub,axis=0)
    return cols


def solve_matrix_inf(A,b,check=False):
    cols = linear_dept_cols(A)
    L,U = linalg.lu(A,permute_l = True)
    U2 = np.transpose(U)[np.where(cols!=0)]
    U2 = np.transpose(U2)
    siz = len(np.where(cols!=0)[0])
    b2 = np.matmul(np.linalg.inv(L),b)
    U3 = U2[:siz]
    b3 = b2[np.where(cols!=0)]
    if np.all(b2[np.where(cols==0)]==0):
        soln1 = np.linalg.solve(U3,b3)
        soln2 = np.zeros(len(b))
        soln2[np.where(cols!=0)] = soln1
        if check:
            if np.allclose(b,np.matmul(A,soln2)):
                print("valid")
        return soln2
    
def solve_matrix(A,b):
    if np.abs(np.linalg.det(A)) > 1e-9:
        return np.linalg.solve(A,b)
    else:
        return solve_matrix_inf(A,b)
