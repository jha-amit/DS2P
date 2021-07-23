

# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 16:46:26 2020

@author: amit
"""
import sys
import numpy as np
import random
import math
import matplotlib as mpl
import numpy as np
from scipy.special import comb
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from pygem import IDW
from statistics import mean
import time
from scipy.integrate import simps
import xlsxwriter

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    

def Patchpoints(nTimes,n,l,theta):
    Patchx=np.ones((n,n))
    Patchy=np.ones((n,n))
    Patchz=np.ones((n,n))
    
    Patchx[:,:]=np.array([np.linspace(j,j+n-1,n) for j in\
                                      range(0,n)])
        
    Patchy[:,:]=np.array([np.linspace(-j,-j+n-1,n) for j in\
                                      range(0,n)])
        
    # This rotation matrix is transpose of the roatation matrix that rotates counterclockwise. 
    # Then my point that i want to roatate will be a row vector or[X,Y] type.
    # Also the angle theta that we are giving input is 
    Rotation_matrix=np.array([[math.cos(theta),math.sin(theta)],[-math.sin(theta),math.cos(theta)]])
    for i in range(0,n):
        for j in range (0,n):
            [Patchx[i,j],Patchy[i,j]]=np.dot([Patchx[i,j],Patchy[i,j]],Rotation_matrix)
    
    # mux=Patchx[2::n,3]
    # muy=Patchy[2::n,3]
    Patchx=Patchx/l
    Patchy=Patchy/l
    # WE CAN PROVIDE A FUNCTION INPUT BY USER TO CALCULATE THE COST.
    # BUT HERE WE ARE PROVIDING COST DIRECTLY FROM USER.
    # mux=Patchx[3,3]
    # muy=Patchy[3,3]   
    # k=0
    # for i in range(0,n):
    #     # if i%3==0 and i>=3:
    #     #         k=k+1
    #     for j in range(0,n):
            
    #         #Patchz[i,j]=math.exp(-(Patchx[i,j]-mux[k])**2-(Patchy[i,j]-muy[k])**2)
    #         Patchz[i,j]=math.exp(-(Patchx[i,j]-mux)**2-(Patchy[i,j]-muy)**2)
            
    return Patchx,Patchy
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# HASH FUNCTION
def hash1(n,K):
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
    return Hash

def interpolation(nTimes,n,l,Patchx,Patchy,Patchz):
    
    def bernstein_poly(i, n, t):
        return comb(n, i) * (t**(n - i)) * (1 - t)**i
    k=0
    nPointsC = 3 #Number of points in each row in a patch
    
    #number of rows and columns in interpolated points array
    if n%2!=0:   
        r=nTimes*(n-1)//2-(n-1)//2+1
    else:
        r= nTimes*((n-1)//2)-((n-1)//2)+1+(nTimes//2)
        #print(r)
        
    #N = (n-1)//2
    # X,Y,Z values are stacked for each 4X4 patch
    
    xvals=np.zeros((r,r)) 
    yvals=np.zeros((r,r)) 
    zvals=np.zeros((r,r))
    vals1=np.zeros((nTimes,nTimes))
 
#parameters for x and y
    t=np.linspace(0.0, 1.0, nTimes)
#print(xPoints)
    polynomial_arrayV =np.array(
    [bernstein_poly(i, nPointsC - 1, t) for i in range(0, nPointsC)])
    j=0    
    while j<=r-nTimes:
        #print(j,k,r)
        
        xvals[0:nTimes,j:j+nTimes]=np.dot(np.transpose(polynomial_arrayV),\
                                    np.dot(Patchx[0:3,k:k+3],polynomial_arrayV))
        
        xvals[0:nTimes,j:j+nTimes]=np.rot90(np.rot90(xvals[0:nTimes,j:j+nTimes],1),1)
        
        xvals[r-nTimes:r,j:j+nTimes]=np.dot(np.transpose(polynomial_arrayV),\
                                    np.dot(Patchx[n-3:n,k:k+3],polynomial_arrayV))
        
       
        xvals[r-nTimes:r,j:j+nTimes]=np.rot90(np.rot90(xvals[r-nTimes:r,j:j+nTimes],1),1)
               
        
        yvals[0:nTimes,j:j+nTimes]=np.dot(np.transpose(polynomial_arrayV),\
                                    np.dot(Patchy[0:3,k:k+3],polynomial_arrayV))
        
        yvals[0:nTimes,j:j+nTimes]=np.rot90(np.rot90(yvals[0:nTimes,j:j+nTimes],1),1) 
       
        yvals[r-nTimes:r,j:j+nTimes]=np.dot(np.transpose(polynomial_arrayV),\
                                    np.dot(Patchy[n-3:n,k:k+3],polynomial_arrayV))
            
        yvals[r-nTimes:r,j:j+nTimes]=np.rot90(np.rot90(yvals[r-nTimes:r,j:j+nTimes],1),1)
            
        zvals[0:nTimes,j:j+nTimes]=np.dot(np.transpose(polynomial_arrayV),\
                                    np.dot(Patchz[0:3,k:k+3],polynomial_arrayV))
        zvals[0:nTimes,j:j+nTimes]=np.rot90(np.rot90(zvals[0:nTimes,j:j+nTimes],1),1)
        
        zvals[r-nTimes:r,j:j+nTimes]=np.dot(np.transpose(polynomial_arrayV),\
                                    np.dot(Patchz[n-3:n,k:k+3],polynomial_arrayV))
        
        zvals[r-nTimes:r,j:j+nTimes]=np.rot90(np.rot90(zvals[r-nTimes:r,j:j+nTimes],1),1)
              
        j=j+nTimes//2
       
        k=k+1            
   
    k=1
    c=0
    j=0
    i=0
    # calculating the interpolation points by moving a 3x3 patch along
    # 3 columns 
    while j <=r-nTimes:
        k=1
        i=nTimes//2
        #print(i,j,k,c)
        while i<=r-nTimes:            
           
            vals1[:,:] = np.dot(np.transpose(polynomial_arrayV),\
                                  np.dot(Patchx[k:k+3,c:c+3],polynomial_arrayV))
            
            vals1[:,:]=np.rot90(np.rot90(vals1[:,:],1),1)
             
            xvals[i:i+nTimes//2+1,j:j+nTimes] = (vals1[0:(nTimes//2)+1,:]\
                                                +xvals[i:i+nTimes//2+1,j:j+nTimes])/2
                
            xvals[i+nTimes//2:i+nTimes,j:j+nTimes] = vals1[nTimes//2:nTimes,:]
            
       
            vals1[:,:] = np.dot(np.transpose(polynomial_arrayV),\
                                  np.dot(Patchy[k:k+3,c:c+3],polynomial_arrayV))
            
            vals1[:,:]=np.rot90(np.rot90(vals1[:,:],1),1)
            
            yvals[i:i+nTimes//2+1,j:j+nTimes] = (vals1[0:(nTimes//2)+1,:]\
                                                +yvals[i:i+nTimes//2+1,j:j+nTimes])/2
                
            yvals[i+nTimes//2:i+nTimes,j:j+nTimes] = vals1[nTimes//2:nTimes,:]
            
            
            vals1[:,:] = np.dot(np.transpose(polynomial_arrayV),\
                                  np.dot(Patchz[k:k+3,c:c+3],polynomial_arrayV))
            
            vals1[:,:]=np.rot90(np.rot90(vals1[:,:],1),1)    
          
            zvals[i:i+nTimes//2+1,j:j+nTimes] = (vals1[0:(nTimes//2)+1,:]\
                                                +zvals[i:i+nTimes//2+1,j:j+nTimes])/2
            
            zvals[i+nTimes//2:i+nTimes,j:j+nTimes] = vals1[nTimes//2:nTimes,:]
            i=i+nTimes//2
            k=k+1        
       
        c=c+1
        j=j+nTimes//2           
        k=1
    return xvals,yvals,zvals
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
# COST FOR SPARSE CASE
def Cost_diamomdgraph(l,n,zvals,nTimes):
    def simpson(z,nTimes,edgelength):
        if len(z)%2==0:
            dx=3*edgelength/(nTimes-2) #step length
            A=z[0]+z[len(z)-2]
            f=0
        # followed wikipedia here for formula.
            for i in range(1, len(z)-2):
                if i%2==0:
                    f=f+dx*(2*z[i])/3
                else:
                    f=f+dx*(4*z[i])/3
            
        else:
            dx=3*edgelength/(nTimes-1) #step length
            A=z[0]+z[len(z)-1]
            f=0
            for i in range(1, len(z)-1):
                if i%2==0:
                    f=f+dx*(2*z[i])/3
                else:
                    f=f+dx*(4*z[i])/3
                    
            #simpson integral
                
        F=A*dx/3+f 
        return F   
    Cost_horizontal=np.zeros((n,n))
    Cost_vertical=np.zeros((n,n))
    Cost_diag=np.zeros((n,n))
    '''if n%2!=0:   
        r=nTimes*(n-1)//2-(n-1)//2+1
    else:
        r= nTimes*((n-1)//2)-((n-1)//2)+1+(nTimes//2)
        #print(r)'''
    
    edgelength=math.sqrt(2)/l #we are taking the edge length as sqrt(2)
                              #eucledean distance between two nodes. 
    k=0
    c=0
    k1=0
    c1=0
    for i in range(0,n):
        for j in range(1,n):
            # WE ARE TAKING nTimes//2 BECAUSE nTimes is 2 nodes span.
            Cost_horizontal[i,j]=simpson(zvals[k,c:c+(nTimes//2)+1],(nTimes//2)+1,edgelength)
            Cost_vertical[j,i]=simpson(zvals[k1:k1+(nTimes//2)+1,c1],(nTimes//2)+1,edgelength)
            c=c+nTimes//2
            k1=k1+nTimes//2
        k1=0
        c=0
        k=k+nTimes//2
        c1=c1+nTimes//2
   # DIAGONALS OF THE FRAME WILL HAVE A LENGTH OF 2 WHILE HAVING SAME NUMBER OF POINTS.    
    k=0
    c=0
    edgelength_DIAGONAL=math.sqrt(4)/l
    for i in range(1,n-1):
        for j in range(1,n-1):
            Cost_diag[i,j]=simpson([zvals[k+i,c+i] for i in range(0,nTimes//2+1)],nTimes//2+1, edgelength_DIAGONAL) 
            c=c+nTimes//2
        k=k+nTimes//2
        c=0     
    return Cost_horizontal,Cost_vertical,Cost_diag
#  COST FOR THE COMPLETE GRAPH NEAR START.
def complete_graph(l,n,zvals,nTimes):
    def simpson(z,nTimes,edgelength):
        if len(z)%2==0:
            dx=3*edgelength/(nTimes-2) #step length
            A=z[0]+z[len(z)-2]
            f=0
        # followed wikipedia here for formula.
            for i in range(1, len(z)-2):
                if i%2==0:
                    f=f+dx*(2*z[i])/3
                else:
                    f=f+dx*(4*z[i])/3
            
        else:
            dx=3*edgelength/(nTimes-1) #step length
            A=z[0]+z[len(z)-1]
            f=0
            for i in range(1, len(z)-1):
                if i%2==0:
                    f=f+dx*(2*z[i])/3
                else:
                    f=f+dx*(4*z[i])/3
                    
            #simpson integral
                
        F=A*dx/3+f 
        return F   
    Cost=np.zeros((n**2,n**2))
    #Cost_right=np.zeros((n**2,n**2))
    c1=0
    for i in range(0,n**2):
        for j in range(i,n**2):
            edgelength=math.sqrt((l*(j//n-i//n))**2+(l*(j%n-i%n))**2)
            if (j//n-i//n)!=0 and (j%n-i%n)!=0:
                
                I=int((j%n-i%n)/(j//n-i//n))
                
                
                
                Cost[i,j]=simpson([zvals[(i//n)*(nTimes//2)+c,(i%n)*(nTimes//2)+I*c] \
                        for c in range(0,((j//n)-(i//n))*(nTimes//2)+1)],\
                        ((j//n)-(i//n))*(nTimes//2)+1, edgelength)
                Cost[j,i]=Cost[i,j]

                # Cost[i,j]=simpson([zvals[c,(i%n)*(nTimes//2)+I*c] \
                #         for c in np.arange((i//n)*(nTimes//2),(j//n)*(nTimes//2)+1,int((j//n-i//n)/(j%n-i%n)))],\
                #         ((j//n)-(i//n))*(nTimes//2)+1, edgelength)
            # elif (j//n-i//n)==0 and (j%n-i%n)==0:
            #     Cost[i,j]=0
            #     Cost[i,j]=0
                
            elif (j//n-i//n)==0:
                c1=(i//n)*(nTimes//2)
                Cost[i,j]=simpson([zvals[c1,c]\
                        for c in range((i%n)*(nTimes//2),(j%n)*(nTimes//2)+1)],\
                            ((j%n)-(i%n))*(nTimes//2+1), edgelength)
                Cost[j,i]=Cost[i,j]
                
                
            elif (j%n-i%n)==0:
                c1=(i%n)*(nTimes//2)
                Cost[i,j]=simpson([zvals[c,c1]\
                        for c in range((i//n)*(nTimes//2),(j//n)*(nTimes//2)+1)],\
                            ((j//n)-(i//n))*(nTimes//2+1), edgelength)
                Cost[j,i]=Cost[i,j]
                
    return Cost  


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def bellman_ford(N,start_node,target_node,B,highvalue):
    node_value = [highvalue]*N
    node_value[start_node]=0    
    #% Run relaxation for N-1 time 
    for n in range(0,N-1):
        for i in range(0,N):
            for j in range(0,N):
                if B[i,j]!=math.inf: #% only performs operation if there is a path cost provided
                   node_value[j]=min(node_value[i]+B[i,j],node_value[j])
                else:
                   node_value[j]=node_value[j]
                
    #target_cost=node_value[0,target_node] #%get the target node path cost
    k=0
    L=[]
    L.append(target_node)
    while (target_node != start_node):
        for k in range(0,N):           
            if node_value[target_node] == node_value[k]+ B[k,target_node] and B[k,target_node]!=0:
                L.append(k)
                target_node=k                
    return L