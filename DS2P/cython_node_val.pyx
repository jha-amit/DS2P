# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 16:24:53 2021

@author: amit
"""


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
    #long long

#to save us some more time as we are sure about indices not going out of
#bounds and there are no negative index to deal with so we can supress these
#two processes 
@cython.boundscheck(False) # Deactivate bounds checking
@cython.wraparound(False) # Deactivate negative indexing.

#Here we are type declaring the Hash and other variables using fused type function
def node_val(my_type[:, ::1] Hash, my_type[:, ::1] Cost_horizontal, my_type[:, ::1] Cost_vertical, my_type[:, ::1] Cost_diag, my_type[::1] N, my_type K, my_type d):
    
    n = Hash.shape[0]
    k = Hash.shape[1]

    #assert tuple(.shape) == tuple(array_2.shape)
#Here accprding to 
    
    #elif my_type is cython.longlong:
        #dtype = np.longlong
    
    if my_type is int:
        dtype = np.intc
    elif my_type is double:
        dtype = np.double
        
    Node_val = np.random.rand(n, k)*np.inf
    
    Node_val[0,0]=0.0
    Node_val[0,1]=0.0
    Node_val[0,2]=0.0
    cdef double[:,::1] node_val1 = Node_val
    # node_val1[0,0]=0.0 #declaring node_val1[0,:]=[0,0,0] gives error ??
    # node_val1[0,1]=0.0
    # node_val1[0,2]=0.0
    
    # if my_type is int:
    #     dtype = np.intc
    # elif my_type is double:
    #     dtype = np.double
    #print(Cost_horizontal[0,0])
    cdef double  ip, jp,i1, i2,j,l1,l2
    cdef double Cijk=0
    cdef int v1, v2,j1,j2,j3,source_neighbour,sink_neighbour
    source_neighbour=int(d/2)
    sink_neighbour=int(K-source_neighbour)
    j1=1    
    l1=0
    l2=0
    Cijk=0    
    for v1 in range(0,n-1):
        j=(Hash[v1,0])
        if j1!=j: # to save the rest of the calculation in the same layer.
            j1=int(j)
            j2=j1+1
            j3=j1+2
            if Hash[v1,0]==0:
                l1=l1+N[j1]            
                l2=l1+N[j2]
            elif Hash[v1,0]<K-1:            
                l1=l1+N[j1]
                l2=l1+N[j2]+N[j3]
            else:
                l1=l1+N[j1]            
                l2=l1+N[j2]        
        
        
    
        for v2 in range(int(l1),int(l2)):
            #print(l1,l2)
            j=Hash[v2,0]
            #is there an edge from v1 to v2, namely
            #determine if an edge exist to j-level vertex v2
             #from the previous(j-1)-level vertex v1
              #at least the vertices are at the right j-levels
                #get the j coordinate of vertex v2
                #potentially there is an edge, but need to inevstigate further
                #namely, check if it is
                #too early/late to differentiate the paths OR
                #the vertices are far apart (i-distance is at least 2)
            if ((j <source_neighbour and (abs(Hash[v2,0]-Hash[v1,0]) <= 1) and (abs(Hash[v2,1]-Hash[v1,1]) <= 1) and\
                (abs(Hash[v2,2]-Hash[v1,2]) <=1)) or (j > sink_neighbour and\
                (abs(Hash[v2,0]-Hash[v1,0]) <= 1) and (abs(Hash[v2,1]-Hash[v1,1]) <= 1) and\
                (abs(Hash[v2,2]-Hash[v1,2]) <=1)) or\
                ((abs(Hash[v2,1]-Hash[v2,2]) >= d) and\
                (abs(Hash[v2,1]-Hash[v1,1]) <= 1) and\
                (abs(Hash[v2,2]-Hash[v1,2]) <=1) and (abs(Hash[v2,0]-Hash[v1,0]) <= 1))):
                #in this case, add the edge
                #determine its (cumulative) cost first
                i1 = (Hash[v2,1])
                i2 = (Hash[v2,2])
                ip = 0.5*(j-i1)
                jp = 0.5*(j+i1)
                
                if Hash[v1,1]<=Hash[v2,1]: 
                    Cijk = Cost_horizontal[int(ip),int(jp)]

                else:
                    Cijk = Cost_vertical[int(ip),int(jp)]
                

                ip = 0.5*(j-i2)
                jp = 0.5*(j+i2)
                #print(ip,jp)
                if Hash[v1,2]<=Hash[v2,2]: 
                    Cijk = Cost_horizontal[int(ip),int(jp)]+Cijk # changing the inquality here
                else:
                    Cijk = Cost_vertical[int(ip),int(jp)]+Cijk
                
                  
                if (node_val1[v2,0]>(Cijk+node_val1[v1,0])):           
                    node_val1[v2,0]=Cijk+node_val1[v1,0]
                    node_val1[v2,1]=v1
                    node_val1[v2,2]=v2

                #node_val1[v2,:]=[(Cijk+node_val1[v1,0]),(v1),(v2)]
    return Node_val

    
    