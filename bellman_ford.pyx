# -*- coding: utf-8 -*-
"""
Created on Mon Sep 14 16:46:18 2020

@author: amit
"""
import sys
import numpy as np
import random
import numpy as np
import random
import math
cimport cython
import time
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
 
ctypedef fused my_type:
    int
    double
    long

#to save us some more time as we are sure about indices not going out of
#bounds and there are no negative index to deal with so we can supress these
#two processes 
@cython.boundscheck(False) # Deactivate bounds checking
@cython.wraparound(False) # Deactivate negative indexing.

#Here we are type declaring the Hash and other variables using fused type function

def bellman_ford(my_type[:, ::1] B, my_type N, my_type start, my_type target, my_type highvalue):

    #assert tuple(.shape) == tuple(array_2.shape)
#Here accprding to 
    
    #elif my_type is cython.longlong:
        #dtype = np.longlong
    if my_type is int:
        dtype = np.intc
    elif my_type is double:
        dtype = np.double
    elif my_type is long:
        dtype = np.double

    
    cdef int n, i, j, k
    
    node_value1= np.ones(int(N))*highvalue
    node_value1[int(start)]=0.0
    cdef double[::1] node_value = node_value1
    
    
    #B[0,1:]=highvalue    
    #% Run relaxation for N-1 time 
    for n in range(0,int(N-1)):
        for i in range(0,int(N)):
            for j in range(0,int(N)):
                if B[i,j]!=np.inf: #% only performs operation if there is a path cost provided
                   node_value[j]=min(node_value[i]+B[i,j],node_value[j])
                   #min(node_value[i]+B[i,j],node_value[j])
                   #print(node_value[j])
                else:
                   node_value[j]=node_value[j]
             
    #target_cost=node_value[0,target_node] #%get the target node path cost
    k=0
    cdef list L= []    
    L.append(target)
    L.append(node_value[int(target)])
    #target_original=target
    while (target != start):
        print(target,k,node_value[int(target)],node_value[k])
        k=0
        for k in range(int(N)):                
            if (node_value[int(target)] == node_value[k]+ B[k,int(target)]) and (B[k,int(target)]!=0.0):
                L.append(k)                
                L.append(node_value[k])                            
                target=k                      
    return L

    
    