# -*- coding: utf-8 -*-
"""
Created on Mon Sep 14 16:59:18 2020

@author: amit
"""
import sys
#cimport cython
import numpy as np
import random
import numpy as np
import math
import time
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
from cython_node_val import node_val

start_time=time.time()
def shortest_path(K):
    if ((K%2!=0) or (K < 0)):
        print("K must be even and positive!")
    else:
        C = np.random.rand(int(K/2+1), int(K/2+1))
    
        #number of nodes 
        n=int(1/48*((K+2)*(K+4)*(2*K+6)))
        #Node co-ordinates 
        Hash=np.zeros((n,3))
    
        #first Nodes first co-ordinate
        Hash[0][:]=[0,0,0]
        j=1
        n=1
        count=0
        #Populate nodes according to the topography
        while n <=1/48*((K+4)*(K+2)*(K+6))-1:
            if count==j+1:
                j=j+1
                count=0
            i2=-j+2*count
            i1=np.linspace(i2,j,j-count+1) 
            for i in range(0,len(i1)):        
                Hash[n+i,0]=j
                Hash[n+i,1]=i1[i]
                Hash[n+i,2]=i2
            count=count+1
            n=n+len(i1)
        count=0
        j=j+1
        while n <=1/48*((K+2)*(K+4)*(2*K+6))-1:
            if count==K-j+1:
                j=j+1
                count=0
            i2=-(K-j)+2*count
            i1=np.linspace(i2,K-j,K-j-count+1)
            for i in range(0,len(i1)):        
                Hash[n+i,0]=j
                Hash[n+i,1]=i1[i]
                Hash[n+i,2] = i2          
            n=n+len(i1)
            count=count+1
            
    Node_val=node_val(Hash,C,K)
    end_time=time.time()
    print(end_time-start_time)
    sPath=np.zeros((K,3))
    sPath[K-1,:]=Node_val[n-1,:]
    for m in range(K-2,-1,-1):
        #print(m)
        sPath[m,:]=Node_val[int(sPath[m+1,1])]
    
    print("Amit")