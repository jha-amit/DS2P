

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

#from test_function1 import Patch_runningaverage,bernstein_poly,simpson,\
#    simpson3,eight_point_rule,COST
#from bezier_quadratic import bezier_curve_quadratic
#from Bezier_movingaveragetest import Bezier_quadratic
#from cython_bezier import bezier_curve_quadratic


#K=30
l=5
nTimes=11
#n=100
#n=int(12/2)+1 #number of rows or column. This is fixed like this because so as to 

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
 #nTimes should be such that nTimes%3==1 "this needs to change"
def bezier_curve_quadratic(nTimes,n,l,theta):
    Patchx=np.ones((n,n))
    Patchy=np.ones((n,n))
    Patchz=np.ones((n,n))
    
    def bernstein_poly(i, n, t):
        return comb(n, i) * (t**(n - i)) * (1 - t)**i

    
    Patchx[:,:]=np.array([np.linspace(j,j+n-1,n) for j in\
                                      range(0,n)])
        
    Patchy[:,:]=np.array([np.linspace(-j,-j+n-1,n) for j in\
                                      range(0,n)])
        
    #Rotation_matrix=np.array([[math.cos()]])
    Rotation_matrix=np.array([[math.cos(theta),-math.sin(theta)],[math.sin(theta),math.cos(theta)]])
    for i in range(0,n):
        for j in range (0,n):
            [Patchx[i,j],Patchy[i,j]]=np.dot([Patchx[i,j],Patchy[i,j]],Rotation_matrix)
    
    # mux=Patchx[2::n,3]
    # muy=Patchy[2::n,3]
    Patchx=Patchx/l
    Patchy=Patchy/l
    
    mux=Patchx[3,3]
    muy=Patchy[3,3]   
    k=0
    for i in range(0,n):
        # if i%3==0 and i>=3:
        #         k=k+1
        for j in range(0,n):
            
            #Patchz[i,j]=math.exp(-(Patchx[i,j]-mux[k])**2-(Patchy[i,j]-muy[k])**2)
            Patchz[i,j]=math.exp(-(Patchx[i,j]-mux)**2-(Patchy[i,j]-muy)**2)
   
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    k=0
    nPointsC = 3 #Number of points in each row in a patch
    #nPointsy = 2 #Number of points in each column in a patch
    #nMatrix = ((n-1)//2)*((n-1)//2) # Number of 4x4patch
    
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
            #print(i,j,k,c)
            #xvals[i:i+nTimes,j:j+nTimes] = np.dot(np.transpose(polynomial_arrayV),\
            #                      np.dot(Patchx[r:r+3,c:c+3],polynomial_arrayV)) 
            
            #print(i,k)
            vals1[:,:] = np.dot(np.transpose(polynomial_arrayV),\
                                  np.dot(Patchx[k:k+3,c:c+3],polynomial_arrayV))
            
            vals1[:,:]=np.rot90(np.rot90(vals1[:,:],1),1)
            #print(vals1)
            #print(polynomial_arrayV.shape,Patchx[k:k+3,c:c+3].shape,k,r)   
            xvals[i:i+nTimes//2+1,j:j+nTimes] = (vals1[0:(nTimes//2)+1,:]\
                                                +xvals[i:i+nTimes//2+1,j:j+nTimes])/2
                
            xvals[i+nTimes//2:i+nTimes,j:j+nTimes] = vals1[nTimes//2:nTimes,:]
            
            #print(xvals[i+nTimes//2:i+nTimes,j:j+nTimes])
            
            #print(i,vals1[nTimes//2:nTimes,:].shape,xvals[i+nTimes//2:i+nTimes,j:j+nTimes].shape)
            vals1[:,:] = np.dot(np.transpose(polynomial_arrayV),\
                                  np.dot(Patchy[k:k+3,c:c+3],polynomial_arrayV))
            
            vals1[:,:]=np.rot90(np.rot90(vals1[:,:],1),1)
            
            yvals[i:i+nTimes//2+1,j:j+nTimes] = (vals1[0:(nTimes//2)+1,:]\
                                                +yvals[i:i+nTimes//2+1,j:j+nTimes])/2
                
            yvals[i+nTimes//2:i+nTimes,j:j+nTimes] = vals1[nTimes//2:nTimes,:]
            
            #print(i,vals1[nTimes//2:nTimes,:].shape,yvals[i+nTimes//2:i+nTimes,j:j+nTimes].shape)
            
            vals1[:,:] = np.dot(np.transpose(polynomial_arrayV),\
                                  np.dot(Patchz[k:k+3,c:c+3],polynomial_arrayV))
            
            vals1[:,:]=np.rot90(np.rot90(vals1[:,:],1),1)    
            #print(vals1)
            zvals[i:i+nTimes//2+1,j:j+nTimes] = (vals1[0:(nTimes//2)+1,:]\
                                                +zvals[i:i+nTimes//2+1,j:j+nTimes])/2
            
            zvals[i+nTimes//2:i+nTimes,j:j+nTimes] = vals1[nTimes//2:nTimes,:]
            i=i+nTimes//2
            k=k+1
            #print(k)
       
        c=c+1
        j=j+nTimes//2           
        k=1
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
    def simpson(z,nTimes,edgelength):
        if len(z)%2==0:
            dx=3*edgelength/(nTimes-2) #step length
            A=z[0]+z[len(z)-2]
            f=0
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
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
    cost=np.ones((n,n))*9999
    Cd=np.ones((n,n))*9999
    '''if n%2!=0:   
        r=nTimes*(n-1)//2-(n-1)//2+1
    else:
        r= nTimes*((n-1)//2)-((n-1)//2)+1+(nTimes//2)
        #print(r)'''
    
    edgelength=math.sqrt(2)/l
    k=0
    c=0
    k1=0
    c1=0
    for i in range(0,n):
        for j in range(1,n):
            cost[i,j]=simpson(zvals[k,c:c+(nTimes//2)+1],(nTimes//2)+1,edgelength)
            cost[j,i]=simpson(zvals[k1:k1+(nTimes//2)+1,c1],(nTimes//2)+1,edgelength)
            c=c+nTimes//2
            k1=k1+nTimes//2
        k1=0
        c=0
        k=k+nTimes//2
        c1=c1+nTimes//2
    return Patchx,Patchy,Patchz,cost
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Patchx,Patchy,Patchz,cost=bezier_curve_quadratic(nTimes,n,l,math.pi/4)
# start_time=time.time()

# l=3
# cost=np.zeros((n,n))

# nTimes=11
# Patchx,Patchy,Patchz,cost= bezier_curve_quadratic(nTimes,n,l)

# end_time=time.time()
# print(end_time-start_time)
#Xnew,Ynew,Znew = rotation(nTimes)
#Patch1x,Patch1y,Patch1z=Patch_runningaverage(4)   
# X = np.array([[3,5,7,9],[3,5,7,9],[3,5,6,9],[3,5,7,9]])
# Y = np.array([[3,3,3,3],[1,1,1,1],[-1,-1,0,-1],[-3,-3,-3,-3]])
# Z=np.array([[Patchz[0,3],Patchz[1,4],Patchz[2,5],Patchz[3,6]],\
#                   [Patchz[1,2],Patchz[2,3],Patchz[3,4],Patchz[4,5]],\
#                       [Patchz[2,1],Patchz[3,2],Patchz[3,3],Patchz[5,4]],\
#                     [Patchz[3,0],Patchz[4,1],Patchz[5,2],Patchz[6,3]]])
    
 
    
# X=np.array([[Patchx[1,3],Patchx[2,4],Patchx[3,5]],[Patchx[2,2],Patchx[3,3],Patchx[4,4]],\
#            [Patchx[3,1],Patchx[4,2],Patchx[5,3]]])
# Y=np.array([[Patchy[1,3],Patchy[2,4],Patchy[3,5]],[Patchy[2,2],Patchy[3,3],Patchy[4,4]],\
#            [Patchy[3,1],Patchy[4,2],Patchy[5,3]]])
# Z=np.array([[Patchz[1,3],Patchz[2,4],Patchz[3,5]],[Patchz[2,2],Patchz[3,3],Patchz[4,4]],\
#            [Patchz[3,1],Patchz[4,2],Patchz[5,3]]])
#     # use bezier_test function here
# #Xnew,Ynew,Znew = Bezier_quadratic(X, Y, Z, 11,3)
# #Xnew,Ynew,Znew = bezier_curve_quadratic(11,3,l) 

# fig = plt.figure(num=1, clear=True)
# ax = fig.gca(projection='3d')
#ax.scatter3D(nodesx, nodesy, nodesz, cmap='green')
#ax.scatter3D(xvals, yvals, zvals)
#ax.scatter3D(xvals1, yvals1, zvals1)
#ax.scatter3D(Patch1x, Patch1y, Patch1z, cmap='red')
#ax.scatter3D(Patchx, Patchy, Patchz, cmap='red')

#ax.scatter3D(Xnew, Ynew, Znew, cmap='orange')
#ax.scatter3D(X, Y, Z, cmap='yellow')  

#cost function



# cost[0,0]=0
# testcost=np.zeros((4,4))
# #testcost=COST(Xnew,Ynew,Znew,nTimes,4)
# edgelength=2*math.sqrt(2)/l
# #cost_test1=simpson([Znew[i,nTimes//2-i] for i in range(0,nTimes//2+1)], nTimes//2+1, edgelength)
# cost_original=cost[2,3]+cost[2,4]
