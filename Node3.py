import numpy as np
import matplotlib.pyplot as plt
from utilities import *


'''
def is_independent(A,v):
    v = np.array(v)
    A2 = np.concatenate((A,np.transpose(v.reshape(1,-1))),axis=1)
    conj = np.transpose(A2)@A2
    if np.abs(np.linalg.det(conj)) != 0:
        if A2.shape[0] == A2.shape[1]:
            if np.abs(np.linalg.det(A2)) == 0:
                return A,False
        return A2,True
    return A,False
'''


def is_independent(A,v):
    v = np.array(v)
    svdsol = svdsolve(A,v)
    diff = np.abs(A@svdsol - v)
    if np.max(diff) < 1e-2:
        return A,False
    else:
        return np.concatenate((A,np.transpose(v.reshape(1,-1))),axis=1),True

def get_node_soln(n,nodes,dim,D,h):
    diffx = n.get_dim('x') - nodes.get_dim('x')
    diffy = n.get_dim('y') - nodes.get_dim('y')
    dist = np.sqrt(diffx**2 + diffy**2)
    arr = np.argsort(dist)
    A = np.array(make_combos([diffx[arr[0]],diffy[arr[0]] ],D+h))
    A = np.transpose(A.reshape(1,-1))
    inds = [arr[0]]
    #print(A.shape)
    for i in arr[1:]:
        v = np.array(make_combos([diffx[arr[i]],diffy[arr[i]] ],D+h))
        A,stat = is_independent(A,v)
        if stat:
            inds.append(arr[i])
        if A.shape[0] == A.shape[1]:
            break
    thelist = make_combos2('xy',D+h)
    b = np.zeros(len(thelist))
    b[thelist.index(dim*D)] = factorial(D)
    soln = svdsolve(A,b)
    #soln = np.linalg.solve(A,b)
    return soln,inds

def test(n,nodes,dim,D,h):
    diffx = n.get_dim('x') - nodes.get_dim('x')
    diffy = n.get_dim('y') - nodes.get_dim('y')
    dist = np.sqrt(diffx**2 + diffy**2)
    arr = np.argsort(dist)
    A = np.array(make_combos([diffx[arr[0]],diffy[arr[0]] ],D+h))
    A = np.transpose(A.reshape(1,-1))
    inds = [arr[0]]
    #print(A.shape)
    #print(n.get_dim('x'),n.get_dim('y'))
    for i in arr[1:]:
        v = np.array(make_combos([diffx[arr[i]],diffy[arr[i]] ],D+h))
        A,stat = is_independent(A,v)
        if stat:
            inds.append(arr[i])
            #print(v)
            #print(diffx[arr[i]],diffy[arr[i]])
            #print(nodes.get_dim('x')[arr[i]],nodes.get_dim('y')[arr[i]])
        if A.shape[0] == A.shape[1]:
            break
    thelist = make_combos2('xy',D+h)
    b = np.zeros(len(thelist))
    b[thelist.index(dim*D)] = factorial(D)
    soln = svdsolve(A,b)
    return soln,inds,A,b

def make_stiffness_5(nodes,D,h,dim,labelname = None):
    n = nodes
    M = n.lennodes
    A = np.zeros((M,M))
    if D == 0:
        if labelname:
            ind = n.labnames[labelname]
            w = np.where(n.labels == ind)[0]
            for i in w:
                A[i][i] = 1.
            return A
        else:
            return np.eye(M)
    if not labelname:
        for i in range(0,n.lennodes):
            #print(i)
            soln,ids = get_node_soln(n.get_node(i),n,dim,D,h)
            A[i][ids] = soln
    else:
        ind = n.labnames[labelname]
        w = np.where(n.labels == ind)[0]
        for i in w:
            ids, soln = diff_gen(i,n,D,h,dim)
            A[i][ids] = soln

    return A