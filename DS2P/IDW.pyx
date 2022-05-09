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
def IDW(my_type[::1] X, my_type[::1] Y, my_type[::1] Z, my_type[::1] X_i, my_type[::1] Y_i, my_type latlimit, my_type longlimit):

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

    cdef int s, p
    cdef double d,u,suminf
    
    lstxyzi_out = np.zeros(X_i.shape[0])
    cdef double[::1] lstxyzi = lstxyzi_out
    for p in range(len(X_i)):
        sumsup=[]
        z_temp=[]
        for s in range(len(X)):
            if abs(X[s]-X_i[p])<=latlimit and abs(Y[s]-Y_i[p])<=longlimit:
                d = math.sqrt(((X[s]-X_i[p])*111.1)**2+((Y[s]-Y_i[p])*111.32*math.cos(X[s]*3.14/180))**2)                
                if d!=0:
                    sumsup.append(1/d)
                    z_temp.append(Z[s])
        suminf = np.sum(sumsup)
        sumsup = np.sum(np.array(sumsup) * np.array(z_temp))
        if suminf!=0:
            u = sumsup/suminf
            lstxyzi[p]=u        
    return(lstxyzi_out)

    
    