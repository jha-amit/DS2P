

# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 16:46:26 2020

@author: amit
"""
from cmath import inf
import sys
from matplotlib import projections
import numpy as np
import random
import math
import matplotlib as mpl
import numpy as np
from scipy.special import comb
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
#from pygem import IDW
from statistics import mean
import time
from scipy.integrate import simps
import xlsxwriter
import random
import csv

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    

def Patchpoints(n,theta,alpha): # number of nodes per row or column, Theta=angle of rotation for alighning the start and end,
                                # alpha=angle of inclination of lines with horizontal
    Patchx=np.ones((n,n))
    Patchy=np.ones((n,n))
    Patchz=np.ones((n,n))
    
    Patchx[:,:]=np.array([np.linspace(math.tan(alpha)*j,math.tan(alpha)*(j+n-1),n) for j in\
                                      range(0,n)])
        
    Patchy[:,:]=np.array([np.linspace(-j,-j+n-1,n) for j in\
                                      range(0,n)])
        
    # This rotation matrix is transpose of the roatation matrix that rotates counterclockwise. 
    # Then my point that i want to roatate will be a row vector or [X,Y] type.
    # Also the angle theta that we are giving input is 
    Rotation_matrix=np.array([[math.cos(theta),math.sin(theta)],[-math.sin(theta),math.cos(theta)]])
    #print(Patchx,Patchy)
    for i in range(0,n):
        for j in range (0,n):
            [Patchx[i,j],Patchy[i,j]]=np.dot([Patchx[i,j],Patchy[i,j]],Rotation_matrix)
    
    # mux=Patchx[2::n,3]
    # muy=Patchy[2::n,3]
    # WE CAN PROVIDE A FUNCTION INPUT BY USER TO CALCULATE THE COST.

    return Patchx,Patchy
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# HASH FUNCTION
def hash1(n,K): # K=number of layers stating from 0.
#Node co-ordinates 
    Hash=np.zeros((n,3))

    #first Nodes first co-ordinate
    Hash[0,:]=[0,0,0]
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

def interpolation(nTimes,n,Patchz): # nTimes = number of interpolted points per row, 
                                        #n=number of nodes per row or column, Patch=costs at nodes  
    
    def bernstein_poly(i, n, t):
        return comb(n, i) * (t**(i)) * (1 - t)**(n-i)
    k=0
    nPointsC = 2 #Number of points in each row in a patch
    
    #number of rows and columns in interpolated points array

    r=(nTimes-1)*(n-1)+1
    zvals=np.zeros((r,r))
 
#parameters for x and y
    t=np.linspace(0.0, 1.0, nTimes)
#print(xPoints)
    polynomial_arrayV =np.array(
    [bernstein_poly(i, nPointsC - 1, t) for i in range(0, nPointsC)])
    c=0
    j=0

    while j <=r-nTimes+1:
        k=0
        i=0
        while i<=r-nTimes+1:
        
           
            zvals[i:i+nTimes,j:j+nTimes] = np.dot(np.transpose(polynomial_arrayV),\
                                  np.dot(Patchz[k:k+2,c:c+2],polynomial_arrayV))
            
            #vals1[:,:]=np.rot90(np.rot90(vals1[:,:],1),1)    
            #print(vals1)
            i=i+nTimes-1
            k=k+1
        c=c+1
        j=j+nTimes-1           
        k=0

    return zvals
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
# COST FOR SPARSE CASE
def Cost_diamondgraph(n,zvals,nTimes,L,L_Diag):
         #n=number of nodes per row, zvals=interpolated costs,
                                                     # L=edge length,L_Diag=length of longer diagonal.
    def simpson(z,nTimes,edgelength):
        dx=edgelength/(nTimes-1)
        
        if (len(z)-1)%3==0:
             #step length
            A=z[0]+z[len(z)-1]
            f=0
        # followed wikipedia here for formula.
            for i in range(1, len(z)-1):
                if i%3==0:
                    f=f+2*z[i]
                else:
                    f=f+3*z[i]
            
        else:
            print('enter number of divisions divisible by 3')
            
        F=(A+f)*3*dx/8
        return F

    Cost_horizontal=np.zeros((n,n))
    Cost_vertical=np.zeros((n,n))
    Cost_diag=np.zeros((n,n))
    '''if n%2!=0:   
        r=nTimes*(n-1)//2-(n-1)//2+1
    else:
        r= nTimes*((n-1)//2)-((n-1)//2)+1+(nTimes//2)
        #print(r)'''
    
    edgelength=L #we are taking the edge length as sqrt(2)
                              #eucledean distance between two nodes. 
    k=0
    c=0
    k1=0
    c1=0
    for i in range(0,n):
        for j in range(1,n):
            # WE ARE TAKING nTimes//2 BECAUSE nTimes is 2 nodes span.
            Cost_horizontal[i,j]=simpson(zvals[k,c:c+nTimes],nTimes,edgelength)            
            c=c+nTimes-1
         
        c=0
        k=k+nTimes-1
       

    for i in range(1,n):
        for j in range(0,n):
            
            Cost_vertical[i,j]=simpson(zvals[k1:k1+nTimes,c1],nTimes,edgelength)
            c1=c1+nTimes-1
        c1=0
        k1=k1+nTimes-1
   # DIAGONALS OF THE FRAME WILL HAVE A LENGTH OF 2 WHILE HAVING SAME NUMBER OF POINTS.    
    k=0
    c=0
    edgelength_DIAGONAL=L_Diag
    for i in range(1,n):
        for j in range(1,n):
            Cost_diag[i,j]=simpson([zvals[k+i,c+i] for i in range(0,nTimes)],nTimes, edgelength_DIAGONAL) 
            c=c+nTimes-1
        k=k+nTimes-1
        c=0
    Cost_horizontal[:,0]=inf
    Cost_vertical[0,:]=inf
    Cost_diag[:,0]=inf
    Cost_diag[0,:]=inf


    return Cost_horizontal,Cost_vertical,Cost_diag
#  COST FOR THE COMPLETE GRAPH NEAR START.(Input arguments l=for scaling the edge size as per requirement,
#  n= for nxn array of nodes, zvals=interpolated cost values,nTimes=number of interpolated points(nTimesxnTimes))
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
    Cost=np.ones((n**2,n**2))*math.inf
    #Cost_right=np.zeros((n**2,n**2))
    c1=0
    for i in range(0,n**2):
        for j in range(i,n**2):
            edgelength=math.sqrt((l*(j//n-i//n))**2+(l*(j%n-i%n))**2)
            if (j//n-i//n)!=0 and (j%n-i%n)!=0:
                
                I=int((j%n-i%n)/(j//n-i//n))
                
                
                
                Cost[i,j]=simpson([zvals[(i//n)*(nTimes-1)+c,(i%n)*(nTimes-1)+I*c] \
                        for c in range(0,((j//n)-(i//n))*(nTimes-1)+1)],\
                        ((j//n)-(i//n))*(nTimes), edgelength)
                Cost[j,i]=Cost[i,j]

# for elements in same row and same column
                
            elif (j//n-i//n)==0 and (j%n-i%n)!=0:
                c1=(i//n)*(nTimes-1)
                
                Cost[i,j]=simpson([zvals[c1,c]\
                        for c in range((i%n)*(nTimes-1),(j%n)*(nTimes-1)+1)],\
                            ((j%n)-(i%n))*(nTimes), edgelength)
                Cost[j,i]=Cost[i,j]

                
                
            elif (j%n-i%n)==0 and (j//n-i//n)!=0:
                c1=(i%n)*(nTimes-1)
                Cost[i,j]=simpson([zvals[c,c1]\
                        for c in range((i//n)*(nTimes-1),(j//n)*(nTimes-1)+1)],\
                            ((j//n)-(i//n))*(nTimes), edgelength)
                Cost[j,i]=Cost[i,j]
                
    return Cost  


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def bellman_ford(B,N,start_node,target_node,highvalue): # B is cost matrix
    
    node_value = [highvalue]*N
    node_value[start_node]=int(0)
    
    #B[0,1:]=highvalue    
    #% Run relaxation for N-1 time 
    for n in range(0,N-1):
        for i in range(0,N):
            for j in range(0,N):
                if B[i,j]!=math.inf: #% only performs operation if there is a path cost provided
                   node_value[j]=min(node_value[i]+B[i,j],node_value[j])
                   #print(node_value[j])
                else:
                   node_value[j]=node_value[j]
    #print(node_value)
                
    #target_cost=node_value[0,target_node] #%get the target node path cost
    k=0
    L=[]    
    L.append(target_node)
    L.append(node_value[target_node])

    while (target_node != start_node):
        #print(1)
        k=0
        for k in range(N):                
            if node_value[target_node] == node_value[k]+ B[k,target_node] and B[k,target_node]!=0:
                L.append(k)                
                L.append(node_value[k])
                            
                target_node=k
                #print(k,B[k,target_node])

                      
    return L
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

def Cost_radial(n_circum,zvals,nTimes,n_radial,radial_nodesx,radial_nodesy):

            #n=number of nodes per row, zvals=interpolated costs,
                                                        # L=edge length,L_Diag=length of longer diagonal.
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
    Cost=np.ones((n_radial*(n_circum),n_radial*(n_circum)))*inf
    #Cost_right=np.zeros((n**2,n**2))
    c1=0
    
#ax.scatter3D(radial_nodes_endx, radial_nodes_endy, Patchz_start, s=.5,c='red',label='End')
    for i in range(0,(n_radial*n_circum)):
        for j in range(i,(n_radial*n_circum)):

            edgelength=math.sqrt((radial_nodesx[j//n_circum,j%n_circum]-radial_nodesx[i//n_circum,i%n_circum])\
                            **2+(radial_nodesy[j//n_circum,j%n_circum]-radial_nodesy[i//n_circum,i%n_circum])\
                            **2)
            
            
            if edgelength==0:
                Cost[i,j]=math.inf

            elif (j//n_circum-i//n_circum)!=0 and (j%n_circum-i%n_circum)!=0:
                
                I=int((j%n_circum-i%n_circum)/(j//n_circum-i//n_circum))               
                #print(edgelength)
                Cost[i,j]=simpson([zvals[(i//n_circum)*(nTimes-1)+c,(i%n_circum)*(nTimes-1)+I*c] \
                        for c in range(0,((j//n_circum)-(i//n_circum))*(nTimes-1)+1)],\
                        ((j//n_circum)-(i//n_circum))*(nTimes-1)+1, edgelength)

                Cost[j,i]=Cost[i,j]
            
                
# for elements in same row and same column
                
            elif (j//n_circum-i//n_circum)==0 and (j%n_circum-i%n_circum)!=0:
                c1=(i//n_circum)*(nTimes-1)
                
                Cost[i,j]=simpson([zvals[c1,c]\
                        for c in range((i%n_circum)*(nTimes-1),(j%n_circum)*(nTimes-1)+1)],\
                            ((j%n_circum)-(i%n_circum))*(nTimes-1)+1, edgelength)


                Cost[j,i]=Cost[i,j]

                
                
            elif (j%n_circum-i%n_circum)==0 and (j//n_circum-i//n_circum)!=0:
                c1=(i%n_circum)*(nTimes-1)
                #print(c1,i,j)
                Cost[i,j]=simpson([zvals[c,c1]\
                        for c in range((i//n_circum)*(nTimes-1),(j//n_circum)*(nTimes-1)+1)],\
                            ((j//n_circum)-(i//n_circum))*(nTimes-1)+1, edgelength)

                

                #print(Cost[i,j])

                Cost[j,i]=Cost[i,j]
                #print(c1,i,j)
                
    return Cost

def interpolation_radial(nTimes,zvals,Patchz): # nTimes = number of interpolted points per row, 
                                        #n=number of nodes per row or column, Patch=costs at nodes  
    
    def bernstein_poly(i, n, t):
        return comb(n, i) * (t**(i)) * (1 - t)**(n-i)
    k=0
    nPointsC = 2 #Number of points in each row in a patch
    
    #number of rows and columns in interpolated points array

#parameters for x and y
    t=np.linspace(0.0, 1.0, nTimes)
#print(xPoints)
    polynomial_arrayV =np.array(
    [bernstein_poly(i, nPointsC - 1, t) for i in range(0, nPointsC)])
    c=0
    j=0
    #print(polynomial_arrayV)
    while j <=zvals.shape[1]-nTimes+1:
        k=0
        i=0
        while i<=zvals.shape[0]-nTimes+1:
        
            #print(i,j)
            zvals[i:i+nTimes,j:j+nTimes] = np.dot(np.transpose(polynomial_arrayV),\
                                np.dot(Patchz[k:k+2,c:c+2],polynomial_arrayV))


            # xvals[i:i+nTimes,j:j+nTimes] = np.dot(np.transpose(polynomial_arrayV),\
            #                     np.dot(Patchx[k:k+2,c:c+2],polynomial_arrayV))

            # yvals[i:i+nTimes,j:j+nTimes] = np.dot(np.transpose(polynomial_arrayV),\
            #                     np.dot(Patchy[k:k+2,c:c+2],polynomial_arrayV))
            
            #vals1[:,:]=np.rot90(np.rot90(vals1[:,:],1),1)    
            #print(vals1)
            i=i+nTimes-1
            k=k+1
        c=c+1
        j=j+nTimes-1           
        k=0

    return zvals

    #++++++++++++++++++++++++++++++++++++++++++++++

def Patch_terminal(center_start_lat,center_start_long,center_end_lat,center_end_long, n_circum_start, n_radial_start,theta,user_defined_radius_start,\
                n_circum_end, n_radial_end,user_defined_radius_end,conversion):

    # first angle starts from the first point of path1 

    # first_angle_start_path1 = theta +  math.asin(delta/(2*user_defined_radius_start))

    # first_angle_start_path2 = theta -  math.asin(delta/(2*user_defined_radius_start))

    # first_angle_end = theta + math.acos(delta/(2*user_defined_radius_end)) # points to the Path1 start_node

    # Path1 of terminal graph

    radial_nodes_starty=np.array([[r*math.sin(gamma)for gamma in\
                        (np.linspace(theta,2*math.pi+theta,n_circum_start))]\
                             for r in (np.linspace(0,user_defined_radius_start,n_radial_start))])

    radial_nodes_startx=np.array([[r*math.cos(gamma)for gamma in \
                            (np.linspace(theta,2*math.pi+theta,n_circum_start))]\
                                 for r in (np.linspace(0,user_defined_radius_start,n_radial_start))])

    radial_nodes_endy=np.array([[r*math.sin(gamma)for gamma in\
                        (np.linspace(math.pi+theta,3*math.pi+theta,n_circum_end))]\
                             for r in (np.linspace(0,user_defined_radius_end,n_radial_end))])

    radial_nodes_endx=np.array([[r*math.cos(gamma)for gamma in \
                            (np.linspace(theta+math.pi,3*math.pi+theta,n_circum_end))]\
                                 for r in (np.linspace(0,user_defined_radius_end,n_radial_end))])

    # Path 2 of terminal graph

    # radial_nodes_startx_path2=np.array([[r*math.cos(gamma)for gamma in \
    #                         (np.linspace(theta,-math.pi+theta,n_circum_start))]\
    #                              for r in (np.linspace(0,user_defined_radius_start,n_radial_start))])

    # radial_nodes_starty_path2=np.array([[r*math.sin(gamma)for gamma in\
    #                     (np.linspace(theta,-math.pi+theta,n_circum_start))]\
    #                          for r in (np.linspace(0,user_defined_radius_start,n_radial_start))])



    # radial_nodes_endy=np.array([[r*math.sin(gamma)for gamma in (np.linspace(first_angle_end,2*math.pi+first_angle_end,n_circum_end))] for r in (np.linspace(0,user_defined_radius_end,n_radial_end))])
    # radial_nodes_endx=np.array([[r*math.cos(gamma)for gamma in (np.linspace(first_angle_end,2*math.pi+first_angle_end,n_circum_end))] for r in (np.linspace(0,user_defined_radius_end,n_radial_end))])

 # lat long of path1
    radial_nodes_lat_start = center_start_lat + (radial_nodes_starty)*conversion
    radial_nodes_long_start = np.zeros((radial_nodes_lat_start.shape[0],radial_nodes_lat_start.shape[1]))

    for i in range(radial_nodes_lat_start.shape[0]):
         for j in range(radial_nodes_lat_start.shape[1]):
            radial_nodes_long_start[i,j]=  center_start_long + (radial_nodes_startx[i,j])/(111.32*math.cos(radial_nodes_lat_start[i,j]*3.14/180))

    radial_nodes_lat_end = center_end_lat + (radial_nodes_endy)*conversion
    radial_nodes_long_end = np.zeros((radial_nodes_lat_end.shape[0],radial_nodes_lat_end.shape[1]))

    for i in range(radial_nodes_lat_end.shape[0]):
         for j in range(radial_nodes_lat_end.shape[1]):
            radial_nodes_long_end[i,j]=  center_end_long + (radial_nodes_endx[i,j])/(111.32*math.cos(radial_nodes_lat_end[i,j]*3.14/180))

# lat-Long of path2 ---------

    # radial_nodes_lat_start_path2 = center_start_lat + (radial_nodes_starty_path2)*conversion
    # radial_nodes_long_start_path2 = np.zeros((radial_nodes_lat_start_path2.shape[0],radial_nodes_lat_start_path2.shape[1]))

    # for i in range(radial_nodes_lat_start_path2.shape[0]):
    #      for j in range(radial_nodes_lat_start_path2.shape[1]):
    #         radial_nodes_long_start_path2[i,j]=  center_start_long + (radial_nodes_startx_path2[i,j])/(111.32*math.cos(radial_nodes_lat_start_path2[i,j]*3.14/180))

# End terminal graphs 
    # radial_nodes_lat_end = center_end_lat + (radial_nodes_endy)*conversion
    # radial_nodes_long_end = np.zeros((radial_nodes_lat_end.shape[0],radial_nodes_lat_end.shape[1]))

    # for i in range(radial_nodes_lat_end.shape[0]):
    #      for j in range(radial_nodes_lat_end.shape[1]):
    #         radial_nodes_long_end [i,j]=  center_end_long + (radial_nodes_endx[i,j])/(111.32*math.cos(radial_nodes_lat_end[i,j]*3.14/180))


    return  radial_nodes_lat_start, radial_nodes_long_start,\
            radial_nodes_starty, radial_nodes_startx,radial_nodes_lat_end, radial_nodes_long_end,\
            radial_nodes_starty, radial_nodes_startx

def latlong_generator(r,center,input_points,file_name):

    radius_meters=r # radius for boundary of data generation
    #radius_degrees = radius_meters*2*math.pi / 111300
        
    random.seed(10)
    w = [radius_meters*(random.random()) for i in range(input_points)] # variable radius 
    t_forlong =  [math.cos(2*math.pi*(random.random())) for i in range(input_points)]
    t_forlat=[math.sin(2*math.pi*(random.random())) for i in range(input_points)]
    
    temp_long_off=[]
    lat_off=[]


    for number1, number2,number3 in zip(w, t_forlong,t_forlat):        
        lat=center[0]+number1*number3/111.1
        lat_off.append(lat)
        temp_long_off.append(center[1]+number1 * number2/(111.32*math.cos(lat*3.14/180)))

       # print(temp_long_off)

    Longitude_input=temp_long_off
    Latitude_input=lat_off
    random.seed(11)
    Cost_input=[random.random() for i in range(input_points)]
    with open(file_name, "w") as csv_file:
        writer = csv.writer(csv_file, delimiter=',',lineterminator='\n')
        writer.writerow(['Longitude','Latitude','Cost'])
        for i in range(0, len(Latitude_input)):
            writer.writerow([Longitude_input[i],Latitude_input[i],Cost_input[i]])