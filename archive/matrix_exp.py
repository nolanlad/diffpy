import numpy as np

def partition_inv(A):
    '''
    Inversion by partition of A.
    To be used for multiprocessing in the future
    Numerical Recipies: p.70
    could be more efficient
    '''
    N = A.shape[0]
    p = N//2
    P = A[:p,:p]
    S = A[p:,p:]
    Q = A[:p,p:]
    R = A[p:,:p]
    Sinv = None; Pt = None
    if p < 1000:
        Sinv = np.linalg.inv(S)
        Pt = np.linalg.inv((P - (Q@Sinv)@R))
    else:
        Sinv = partition_inv(S)
        Pt = partition_inv((P - (Q@Sinv)@R))
    del S
    del P
    QSinv = Q@Sinv
    del Q
    SinvR = Sinv@R
    del R
    Qt = -1*Pt@(QSinv)
    St = Sinv + SinvR@Pt@QSinv
    del QSinv
    Rt = -1*(SinvR)@Pt
    A[:p,:p] = Pt
    A[p:,p:] = St
    A[:p,p:] = Qt
    A[p:,:p] = Rt
    return A
        
        

A = np.random.random((7001,7001))
Ainv = np.linalg.inv(A)
#Ainvp = partition_inv(A)
#print(np.allclose(Ainv,Ainvp))
