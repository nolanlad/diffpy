import numpy as np 
import matplotlib.pyplot as plt 
from utilities import *
from Node2 import *
from Node3 import *
from pde import Diff,OpList

'''
Solve for steady state solution for heat conduction in a
square of material with adiabatic walls and constant temperatures
at top and bottom
'''
N=25
x = np.linspace(0,1,N)

n = Nodes({"x":x})
labels = np.zeros(n.lennodes)

'''
front = (x == 0.0)
back = (x == 1.0) 
labels[front] = 1; labels[back] = 2

labnam = {"mass":0,"front":1,"back":2}
n.add_labels(labels,labnam)
h = 3
ops = []
ops.append(Diff(n,'mass','x',2,0.0))
ops.append(Diff(n,'back','x',0,0.0))
ops.append(Diff(n,'front','x',2,1.0))

opl = OpList(ops)

opl.combine_linops()

A = opl.lop.mat
b = opl.lop.targ

'''

