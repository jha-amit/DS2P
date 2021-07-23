from django.shortcuts import render
import sys
import numpy as np
import random
import math
import gmplot
import csv
from cython_node_val import node_val
from Cost_modified import interpolation,Cost_diamomdgraph,Patchpoints,hash1,complete_graph,bellman_ford 
import pickle
from .forms import UploadFileForm
from django.conf import settings
from django.core.files.storage import FileSystemStorage
from django.core.files.storage import FileSystemStorage
from django.db import models
from numpy import genfromtxt



def home(request):
    return render(request,'Home.html',{"name":"user"})

def shortest_path1(request):

    conversion=1/111 #conversion from change in latitude to km.
                    #conversion for km to latitude = (111.32*math.cos(lat[k]*math.pi/180)) 
    start_lat=float(request.POST['Start point latitude']) 
    start_long=float(request.POST['Start point longitude'])
    end_lat=float(request.POST['End point latitude'])  
    end_long=float(request.POST['End point longitude'])   
    start=[start_lat,start_long]
    colors = ['red','blue','green','purple','orange','yellow','pink','white']
    gmap5 = gmplot.GoogleMapPlotter(start[0],start[1],13) #assign a map from Gmapplotter 
    gmap5.apikey = request.POST['enter your gmap API key']  
    
    K=int(request.POST['number of layers'])
    if ((K%2!=0) or (K < 0)):
        return render(request,"shortest_path1.html",{'shortest_path1':"K must be an even integer"})
    else:
        l=5     # This parameter is to control the actual edge length. it divides 1 unit in l parts so we
                #can shrink or expand the grid for our uses. here i am using it for optimal size of grid to acheiv
                #  fairly accurate cost interpolations. 
        N=int(K/2)+1
        Patchx=np.zeros((N,N))
        Patchy=np.zeros((N,N))
        Patchz=np.zeros((N,N))  
        nTimes=11
        
        
        # 111 km per degree latitude change and 111.17km for per degree longitude change 
        #theta=0
        theta=math.atan((((end_lat-start_lat)*111)/((end_long-start_long)*111.2))) # need correction for eliptical co-ordinates
        # We need to provide user data as (Lat,Long,Cost) This user data will be in the form of scattered data
        # point. But we need the costs along the edges of diamond frame. I will interpolate 
        # 'Inverse distance weight' method the user input for getting costs at control points which is nodes
        #  of diamond graph.
        print(theta)
        

        Patchx,Patchy= Patchpoints(nTimes,N,l,theta) #cost array
        k=0
        Lat=[1]*N**2
        Long=[1]*N**2
        for i in range(0,N):
            for j in range(0,N):
                Lat[k]=Patchy[i,j]*conversion+start[0]       
                Long[k]=Patchx[i,j]/(111.32*math.cos(Lat[k]*3.14/180))+start[1]
                k=k+1
   
     
# If the user provides (Latitude,Longitude,cost) then we resort to this procedure.
# We will request the post method to provide the labeled data from input.
# Then we convert our patch control points in Latitude[] Longitude[].
# Then we apply IDW algo. to find costs at these latitude longitude.
# give this cost as control point to Bezier curves.
# Get the cost inputs required along edges.
# We can also plot a Cost Heatmap.

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'''
# Using filestorage API provided by Django upload and save 
        Longitude_input=[]
        Latitude_input=[]
        Cost_input=[]
        if request.method == 'POST':
            Cost_matrix = request.FILES['Cost_matrix']
            fs = FileSystemStorage()
            filename = fs.save(Cost_matrix.name, Cost_matrix)

        # with open('C:/Users/amit/Desktop/Thesis codes compilation/Application Python/application_shortest_path/shortest_path/media/Cost_matrix.csv', 'r', encoding='utf-8-sig') as f: 
        #     C = np.genfromtxt(f, dtype=float, delimiter=',')

        # Reading the file from filestorage 
        with open('E:/Project/Thesis_related_codes/Application_Python/application_shortest_path/shortest_path/media/Cost_matrix.csv', newline='') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                Latitude_input.append(float(row['Latitude']))
                Longitude_input.append(float(row['Longitude']))
                Cost_input.append(float(row['Cost']))
                
                        
        # Inverse distance weight algorithm. We get a list of interpolated points for the unknown xi, yi,
        # using known points x, y, z. 
        def idwr(x, y, z, xi, yi):
            lstxyzi = []
            for p in range(len(xi)):
                lstdist = []
                for s in range(len(x)):
                    d = math.sqrt((x[s]-xi[p])**2+(y[s]-yi[p])**2)
                    if d!=0:
                        lstdist.append(d)
                    else:
                        lstdist.append(.001)
                    
                sumsup = list((1 / np.power(lstdist, 2)))
                suminf = np.sum(sumsup)
                sumsup = np.sum(np.array(sumsup) * np.array(z))
                u = sumsup / suminf
                lstxyzi.append(u)        
            return(lstxyzi)
            
        Cost=idwr(Latitude_input,Longitude_input,Cost_input,Lat,Long)

# Write the diamond node costs to a .csv file to enable user for editing  
        with open('E:/Project/Thesis_related_codes/Application_Python/application_shortest_path/shortest_path/media/diamond_frame_nodes.csv', "w") as csv_file:
            writer = csv.writer(csv_file, delimiter=',',lineterminator='\n')
            writer.writerow(['Latitude','Longitude','Cost'])
            for i in range(0, len(Lat)):
                writer.writerow([Lat[i],Long[i],Cost[i]])

    # read it back to the list Cost.

    # By doing this we can make the input interactive.
        Lat=[]
        Long=[]
        Cost=[]
        with open('E:/Project/Thesis_related_codes/Application_Python/application_shortest_path/shortest_path/media/diamond_frame_nodes.csv', newline='') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                Lat.append(float(row['Latitude']))
                Long.append(float(row['Longitude']))
                Cost.append(float(row['Cost']))
        k=0
        # Cost should go as points to be interpolated by bazier curve cost interpolation function.
     
        for i in range(0,N):
            for j in range(0,N):
                Patchz[i,j]=Cost[k]
                k=k+1        
        if N%2!=0:   
            r=nTimes*(N-1)//2-(N-1)//2+1
        else:
            r= nTimes*((N-1)//2)-((N-1)//2)+1+(nTimes//2)
        xvals=np.zeros((r,r)) 
        yvals=np.zeros((r,r)) 
        zvals=np.zeros((r,r))
        xvals,yvals,zvals=interpolation(nTimes,N,l,Patchx,Patchy,Patchz)
        Cost_horizontal,Cost_vertical,Cost_diag=Cost_diamomdgraph(l,N,zvals,nTimes)

    # HERE WE ARE TRYING TO PROVIDE A GRAPH FOR THE START NODE TO MITIGATE THE RESTRICTION NEAR SORCE NODE.
    # test the rotation and scaling of patch
       
    nrect=int(request.POST['graph_at_source']) # number of nodes in the graph near source.
   # GRAPH ON RIGHT HAND SIDE
    Patch1x=np.zeros((nrect,nrect))
    Patch1y=np.zeros((nrect,nrect))
    Patch_leftx=np.zeros((nrect,nrect))
    Patch_lefty=np.zeros((nrect,nrect))
    Patch_leftz=np.zeros((nrect,nrect))
    Patch_rightx=np.zeros((nrect,nrect))
    Patch_righty=np.zeros((nrect,nrect))
    Patch_rightz=np.zeros((nrect,nrect))
    theta1=theta+math.pi/4  
    theta2=theta-3*math.pi/4  
    Rotation_matrix1=np.array([[math.cos(theta1),math.sin(theta1)],[-math.sin(theta1),math.cos(theta1)]])
    Rotation_matrix2=np.array([[math.cos(theta2),math.sin(theta2)],[-math.sin(theta2),math.cos(theta2)]])

    Patch1x[:,:]=np.array([np.arange(0,nrect,1) for i in range(0,nrect)])
    Patch1y[:,:]=np.array([[j for i in range(0,nrect)] for j in range(0,nrect)])
    l1=4*l # this is not the length, This is number of parts a unit length is divided into
    Patch1x=Patch1x/l1
    Patch1y=Patch1y/l1
    for i in range(0,nrect):
        for j in range(0,nrect):
            [Patch_leftx[i,j],Patch_lefty[i,j]]=np.dot([Patch1x[i,j],Patch1y[i,j]],Rotation_matrix1) # this is a clocwise rotation
            [Patch_rightx[i,j],Patch_righty[i,j]]=np.dot([Patch1x[i,j],Patch1y[i,j]],Rotation_matrix2)
    # GRAPH ON LEFT HAND SIDE **** 

    Patch_leftx=Patch_leftx+Patchx[0,0]-math.cos(math.pi/4+theta)/l
    Patch_lefty=Patch_lefty+Patchy[0,0]-math.sin(math.pi/4+theta)/l

    # GRAPH ON RIGHT HAND SIDE  **** 
    
    Patch_rightx=Patch_rightx+Patchx[0,0]-math.cos(math.pi/4-theta)/l
    Patch_righty=Patch_righty+Patchy[0,0]+math.sin(math.pi/4-theta)/l
    #  HERE WE ARE PROVIDING THIS SMALL GRAPTH THE INPUTS REQUIRED FOR CALCULATING SHORTEST PATHS NEAR SOURCE.
    #*** LEFTSIDE GRAPH
    k=0
    Lat_left=[1]*nrect**2
    Long_left=[1]*nrect**2
    Lat_right=[1]*nrect**2
    Long_right=[1]*nrect**2
    for i in range(0,nrect):
        for j in range(0,nrect):
            Lat_left[k]=Patch_lefty[i,j]*conversion+start[0]       
            Long_left[k]=Patch_leftx[i,j]/(111.32*math.cos(Lat_left[k]*3.14/180))+start[1]
            Lat_right[k]=Patch_righty[i,j]*conversion+start[0]       
            Long_right[k]=Patch_rightx[i,j]/(111.32*math.cos(Lat_right[k]*3.14/180))+start[1]
            k=k+1
    # with open('E:/Project/Thesis_related_codes/Application_Python/application_shortest_path/shortest_path/media/left_graph.csv', "w") as csv_file:
    #     writer = csv.writer(csv_file, delimiter=',',lineterminator='\n')
    #     writer.writerow(['Latitude_left','Longitude_left','Cost_left','Latitude_right','Longitude_right','Cost_right'])
    #     for i in range(0, len(Lat_left)):
    #         writer.writerow([Lat_left[i],Long_left[i],' ',Lat_right[i],Long_right[i],' '])

    Longitude_left=[]
    Latitude_left=[]
    Cost_left=[]
    Longitude_right=[]
    Latitude_right=[]
    Cost_right=[]
    if request.method == 'POST':
        Cost_ends = request.FILES['Endgraphs_Cost']
        fs = FileSystemStorage()
        filename = fs.save(Cost_ends.name, Cost_ends)

    with open('E:/Project/Thesis_related_codes/Application_Python/application_shortest_path/shortest_path/media/Endgraphs_Cost_9nodes.csv', newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            Latitude_left.append(float(row['Latitude_left']))
            Longitude_left.append(float(row['Longitude_left']))
            Cost_left.append(float(row['Cost_left']))
            Latitude_right.append(float(row['Latitude_right']))
            Longitude_right.append(float(row['Longitude_right']))
            Cost_right.append(float(row['Cost_right']))

    # INTERPOLATING COSTS FROM INPUT POINTS

    Cost_left=idwr(Latitude_left,Longitude_left,Cost_left,Lat_left,Long_left)
    Cost_right=idwr(Latitude_right,Longitude_right,Cost_right,Lat_right,Long_right)

    # THESE INTERPOLTED COSTS WILL GO INTO THE FINAL INTERPOLATION. HERE THERE WILL BE A
    # DIFFERENT INTERPOLATION FUNCTION IS USED.
    k=0 
    for i in range(0,nrect):
        for j in range(0,nrect):
            Patch_leftz[i,j]=Cost_left[k]
            Patch_rightz[i,j]=Cost_right[k]
            k=k+1    
         
    if nrect%2!=0:   
        r=nTimes*(nrect-1)//2-(nrect-1)//2+1  
    else:
        r= nTimes*((nrect-1)//2)-((nrect-1)//2)+1+(nTimes//2)
    xvals_left=np.zeros((r,r)) 
    yvals_left=np.zeros((r,r)) 
    zvals_left=np.zeros((r,r))
    xvals_right=np.zeros((r,r)) 
    yvals_right=np.zeros((r,r)) 
    zvals_right=np.zeros((r,r))
    Costmatrix_left=np.zeros((nrect**2,nrect**2)) 
    Costmatrix_right=np.zeros((nrect**2,nrect**2))   
    xvals_left,yvals_left,zvals_left=interpolation(nTimes,nrect,l1,Patch_leftx,Patch_lefty,Patch_leftz)
    xvals_right,yvals_right,zvals_right=interpolation(nTimes,nrect,l1,Patch_rightx,Patch_righty,Patch_rightz)
    # THIS GRAPH IS A COMPLETE GRAPH HENCE THE COST MATRIX IS NOT A SPARSE MATRIX
    Costmatrix_left= complete_graph(l1,nrect,zvals_left,nTimes)
    Costmatrix_right= complete_graph(l1,nrect,zvals_right,nTimes)
    
    # WE CALCULATE THE SHORTEST PATH NOW USING BELLMAN-FORD ALGORITHM! 

    #%import cost matrix
    start_row=0
    start_col=nrect//2
    end_row=0
    end_col=nrect-1
    start_node=start_row*nrect+start_col
    target_node=end_row*nrect+end_col
    start_row_right=nrect//2
    start_col_right=0
    end_row_right=nrect-1
    end_col_right=0
    start_node_right=start_row_right*nrect+start_col_right
    target_node_right=end_row_right*nrect+end_col_right
    # B=Costmatrix_left
    # C=Costmatrix_right
    sPath_left=bellman_ford(nrect**2,start_node,target_node,Costmatrix_left,10000000)
    sPath_right=bellman_ford(nrect**2,start_node_right,target_node_right,Costmatrix_right,1000000)

   
    # THIS PART IS FOR THE PURPOSE OF SMALL CORRECTION WE NEED IN THE COST MATRIX.
    # THIS IS TAKING VALUES FROM FRONT END HTML PAGE 

    A=int(request.POST['edit_row'])
    B=int(request.POST['edit_column'])
    D=float(request.POST['cost'])
    #print(C)
    Cost_horizontal[A-1,(B-1)]=D 

    #number of nodes

    n=int(1/48*((K+2)*(K+4)*(2*K+6)))
    Hash=np.zeros((n,3))
    #Node co-ordinates 
    Hash=hash1(n,K)

    N_K=np.array([1.0]*(K+1)) # number of nodes in each layer Note that there are K layers
                                                        #between K+1 planes.
    i=1
    #Number of nodes in a layer
    while i<=K+1:
        if i<=K/2:
            N_K[i-1]=int(i*(i+1)/2)
        else:
            N_K[i-1]=int((K-i+2)*(K-i+3)/2)
        i=i+1

    Node_val=node_val(Hash,Cost_horizontal,Cost_vertical,Cost_diag,N_K,K)

    #print(end_time-start_time)
    sPath=np.zeros((K,3))
    sPath[K-1,:]=Node_val[n-1,:]
    for m in range(K-2,-1,-1):
        sPath[m,:]=Node_val[int(sPath[m+1,1]),:]   
    #print(sPath)
    
    # conversion of 2D node arrangement, Patch ---> list of [latitude, longitude]. 
    for i in range(0,N-1):
        for j in range(0,N-1):           
            lat1=Patchy[i,j]*conversion+start[0]   
            long1=Patchx[i,j]/(111.32*math.cos(lat1*3.14/180))+start[1]
            lat2=Patchy[i,j+1]*conversion+start[0]
            long2=Patchx[i,j+1]/(111.32*math.cos(lat2*3.14/180))+start[1]
            gmap5.plot([lat1,lat2],[long1,long2], colors[2], edge_width = 1.5)
                       
            lat2=Patchy[i+1,j]*conversion+start[0] 
            long2=Patchx[i+1,j]/(111.32*math.cos(lat2*3.14/180))+start[1]
                   
            gmap5.plot([lat1,lat2], [long1,long2], colors[3], edge_width = 1.5)
   
    ends=N-1
    for i in range(0,N-1):
        lat1=Patchy[i,ends]*conversion+start[0]   
        long1=Patchx[i,ends]/(111.32*math.cos(lat1*3.14/180))+start[1]
        lat2=Patchy[i+1,ends]*conversion+start[0]   
        long2=Patchx[i+1,ends]/(111.32*math.cos(lat2*3.14/180))+start[1]
        gmap5.plot([lat1,lat2], [long1,long2], colors[2], edge_width = 1.5)

        lat1=Patchy[ends,i]*conversion+start[0]   
        long1=Patchx[ends,i]/(111.32*math.cos(lat1*3.14/180))+start[1]
        lat2=Patchy[ends,i+1]*conversion+start[0]   
        long2=Patchx[ends,i+1]/(111.32*math.cos(lat2*3.14/180))+start[1]
        gmap5.plot([lat1,lat2], [long1,long2], colors[3], edge_width = 1.5)    
    
    for i in range(0,nrect-1):
        for j in range(0,nrect-1):           
            lat1=Patch_righty[i,j]*conversion+start[0]   
            long1=Patch_rightx[i,j]/(111.32*math.cos(lat1*3.14/180))+start[1]
            lat2=Patch_righty[i,j+1]*conversion+start[0]
            long2=Patch_rightx[i,j+1]/(111.32*math.cos(lat2*3.14/180))+start[1]
            gmap5.plot([lat1,lat2],[long1,long2], colors[2], edge_width = 1.5)
                       
            lat2=Patch_righty[i+1,j]*conversion+start[0] 
            long2=Patch_rightx[i+1,j]/(111.32*math.cos(lat2*3.14/180))+start[1]                   
            gmap5.plot([lat1,lat2], [long1,long2], colors[3], edge_width = 1.5)
   
    ends=nrect-1
    for i in range(0,nrect-1):
        lat1=Patch_righty[i,ends]*conversion+start[0]   
        long1=Patch_rightx[i,ends]/(111.32*math.cos(lat1*3.14/180))+start[1]
        lat2=Patch_righty[i+1,ends]*conversion+start[0]   
        long2=Patch_rightx[i+1,ends]/(111.32*math.cos(lat2*3.14/180))+start[1]
        gmap5.plot([lat1,lat2], [long1,long2], colors[2], edge_width = 1.5)

        lat1=Patch_righty[ends,i]*conversion+start[0]   
        long1=Patch_rightx[ends,i]/(111.32*math.cos(lat1*3.14/180))+start[1]
        lat2=Patch_righty[ends,i+1]*conversion+start[0]   
        long2=Patch_rightx[ends,i+1]/(111.32*math.cos(lat2*3.14/180))+start[1]
        gmap5.plot([lat1,lat2], [long1,long2], colors[3], edge_width = 1.5)

    for i in range(0,nrect-1):
        for j in range(0,nrect-1):           
            lat1=Patch_lefty[i,j]*conversion+start[0]   
            long1=Patch_leftx[i,j]/(111.32*math.cos(lat1*3.14/180))+start[1]
            lat2=Patch_lefty[i,j+1]*conversion+start[0]
            long2=Patch_leftx[i,j+1]/(111.32*math.cos(lat2*3.14/180))+start[1]
            gmap5.plot([lat1,lat2],[long1,long2], colors[2], edge_width = 1.5)
                       
            lat2=Patch_lefty[i+1,j]*conversion+start[0] 
            long2=Patch_leftx[i+1,j]/(111.32*math.cos(lat2*3.14/180))+start[1]                   
            gmap5.plot([lat1,lat2], [long1,long2], colors[3], edge_width = 1.5)
   
    ends=nrect-1
    for i in range(0,nrect-1):
        lat1=Patch_lefty[i,ends]*conversion+start[0]   
        long1=Patch_leftx[i,ends]/(111.32*math.cos(lat1*3.14/180))+start[1]
        lat2=Patch_lefty[i+1,ends]*conversion+start[0]   
        long2=Patch_leftx[i+1,ends]/(111.32*math.cos(lat2*3.14/180))+start[1]
        gmap5.plot([lat1,lat2], [long1,long2], colors[2], edge_width = 1.5)

        lat1=Patch_lefty[ends,i]*conversion+start[0]   
        long1=Patch_leftx[ends,i]/(111.32*math.cos(lat1*3.14/180))+start[1]
        lat2=Patch_lefty[ends,i+1]*conversion+start[0]   
        long2=Patch_leftx[ends,i+1]/(111.32*math.cos(lat2*3.14/180))+start[1]
        gmap5.plot([lat1,lat2], [long1,long2], colors[3], edge_width = 1.5)
    #plotting the shortest route on 2D
    
    # for v1 in range(0,K):        
      
    #     if (sPath[v1,0]>0.0):
    #         [N1,N2] = [int(sPath[v1,1]),int(sPath[v1,2])]
    #         [j1,i11,i12]=[Hash[N1,0],Hash[N1,1],Hash[N1,2]]
    #         [j2,i21,i22]=[Hash[N2,0],Hash[N2,1],Hash[N2,2]]
    #         [X1,Y1,X2,Y2,X3,Y3,X4,Y4]=[int((j1-i11)/2),int((j1+i11)/2),int((j1-i12)/2),int((j1+i12)/2),int((j2-i21)/2),int((j2+i21)/2),int((j2-i22)/2),int((j2+i22)/2)]
            
    #         lat1 = Patchy[X1,Y1]*conversion+start[0]
    #         lat2 = Patchy[X3,Y3]*conversion+start[0]            
    #         long1 = Patchx[X1,Y1]/(111.32*math.cos(lat1*3.14/180))+start[1]
    #         long2 = Patchx[X3,Y3]/(111.32*math.cos(lat2*3.14/180))+start[1]
    #         gmap5.plot([lat1, lat2],[long1, long2],colors[2], edge_width = 2)

    #         lat1 = Patchy[X2,Y2]*conversion+start[0]
    #         lat2 = Patchy[X4,Y4]*conversion+start[0]            
    #         long1 = Patchx[X2,Y2]/(111.32*math.cos(lat1*3.14/180))+start[1]
    #         long2 = Patchx[X4,Y4]/(111.32*math.cos(lat2*3.14/180))+start[1]
    #         gmap5.plot([lat1, lat2],[long1, long2],colors[3], edge_width = 2)
    # k=0
    # for v1 in range(0,len(sPath_left)-1):   
    #     [N1,N2] = [int(sPath_left[v1]),int(sPath_left[v1+1])]
    #     X1=N1//nrect
    #     X2=N2//nrect
    #     Y1=N1%nrect
    #     Y2=N2%nrect       
    #     lat1 = Patch_lefty[X1,Y1]*conversion+start[0]
    #     lat2 = Patch_lefty[X2,Y2]*conversion+start[0]            
    #     long1 = Patch_leftx[X1,Y1]/(111.32*math.cos(lat1*3.14/180))+start[1]
    #     long2 = Patch_leftx[X2,Y2]/(111.32*math.cos(lat2*3.14/180))+start[1]
    #     gmap5.plot([lat1, lat2],[long1, long2],colors[2], edge_width = 2)
    # k=0
    # for v1 in range(0,len(sPath_right)-1):
    #     [N1,N2] = [int(sPath_right[v1]),int(sPath_right[v1+1])]
    #     X1=N1//nrect
    #     X2=N2//nrect
    #     Y1=N1%nrect
    #     Y2=N2%nrect       
    #     lat1 = Patch_righty[X1,Y1]*conversion+start[0]
    #     lat2 = Patch_righty[X2,Y2]*conversion+start[0]            
    #     long1 = Patch_rightx[X1,Y1]/(111.32*math.cos(lat1*3.14/180))+start[1]
    #     long2 = Patch_rightx[X2,Y2]/(111.32*math.cos(lat2*3.14/180))+start[1]
    #     gmap5.plot([lat1, lat2],[long1, long2],colors[3], edge_width = 2)

            
    gmap5.draw( "E:/Project/Thesis_related_codes/Application_Python/application_shortest_path/shortest_path/Template/map.html" )
    # for triggering two different submit buttons
    if request.POST:
        if 'shortest path' in request.POST:
            return render(request,"shortest_path1.html",{'shortest_path1':sPath,'sPath_left':sPath_left})
        elif 'GoogleMaps' in request.POST:
            return render(request,"map.html")
