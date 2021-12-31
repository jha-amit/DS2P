from django.shortcuts import redirect, render
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
import json 
from django.core.files.storage import default_storage

val=None
#Inverse distance interpolation function
def idwr(x, y, z, xi, yi):
    lstxyzi = []
    for p in range(len(xi)):
        lstdist = []
        for s in range(len(x)):
            d = math.sqrt(((x[s]-xi[p])*111.1)**2+((y[s]-yi[p])*111.32)**2)
            if d!=0 and d<=1:
                lstdist.append(d)
            else:
                lstdist.append(100000)
            
        sumsup = list((1 / np.power(lstdist, 2)))
        suminf = np.sum(sumsup)
        sumsup = np.sum(np.array(sumsup) * np.array(z))
        u = sumsup / suminf
        lstxyzi.append(u)        
    return(lstxyzi)

def home(request):
    return render(request,'Home.html',{"name":"user"})



def shortest_path1(request):

    conversion=1/111 #conversion from change in latitude to km.
                    #conversion for km to latitude = (111.32*math.cos(lat[k]*math.pi/180)) 
    start_lat=float(request.POST['Start point latitude']) 
    start_long=float(request.POST['Start point longitude'])
    end_lat=float(request.POST['End point latitude'])  
    end_long=float(request.POST['End point longitude'])
    end=[end_lat,end_long]   
    start=[start_lat,start_long]
    K=int(request.POST['number of layers'])
    nTimes=int(request.POST['Grid density'])
    alpha=float(request.POST['edge_slope'])

    if ((K%2!=0) or (K < 0)):
        return render(request,"shortest_path1.html",{'shortest_path1':"Enter Correct Inputs and Press Input button"})
    else:
        # This parameter is to control the actual edge length. it divides 1 unit in l parts so we
                #can shrink or expand the grid for our uses. here i am using it for optimal size of grid to acheiv
                #  fairly accurate cost interpolations. 
        N=int(K/2)+1
        Patchx=np.zeros((N,N))
        Patchy=np.zeros((N,N))
        Patchz=np.zeros((N,N))
        
        
        # 111 km per degree latitude change and 111.17km for per degree longitude change 
        #theta=0
        theta=math.atan((((end_lat-start_lat)*111)/((end_long-start_long)*111.2))) # need correction for eliptical co-ordinates
        # We need to provide user data as (Lat,Long,Cost) This user data will be in the form of scattered data
        # point. But we need the costs along the edges of diamond frame. I will interpolate 
        # 'Inverse distance weight' method the user input for getting costs at control points which is nodes
        #  of diamond graph.
        #print(theta)
        Patchx,Patchy= Patchpoints(N,theta,alpha) #cost array
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
            myfile = request.FILES['Cost_matrix']
            fs = FileSystemStorage()
            Cost_matrix = fs.save(myfile.name, myfile)

# Using default file storage
#     file = request.FILES['myfile']
#     file_name = default_storage.save(file.name, file)

# #  Reading file from storage
#     file = default_storage.open(file_name)
#     file_url = default_storage.url(file_name)

        # with open('C:/Users/amit/Desktop/Thesis codes compilation/Application Python/application_shortest_path/shortest_path/media/Cost_matrix.csv', 'r', encoding='utf-8-sig') as f: 
        #     C = np.genfromtxt(f, dtype=float, delimiter=',')

        # Reading the file from filestorage 
        with open('C:/Users/Amit/Desktop/Working/application_shortest_path/shortest_path/media/Cost_matrix.csv', newline='') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                Latitude_input.append(float(row['Latitude']))
                Longitude_input.append(float(row['Longitude']))
                Cost_input.append(float(row['Cost']))
                
                        
        # Inverse distance weight algorithm. We get a list of interpolated points for the unknown xi, yi,
        # using known points x, y, z. 
        
            
        Cost=idwr(Latitude_input,Longitude_input,Cost_input,Lat,Long)

# Write the diamond node costs to a .csv file to enable user for editing  
        with open('media/diamond_frame_nodes.csv', "w") as csv_file:
            writer = csv.writer(csv_file, delimiter=',',lineterminator='\n')
            writer.writerow(['Latitude','Longitude','Cost'])
            for i in range(0, len(Lat)):
                writer.writerow([Lat[i],Long[i],Cost[i]])

    # read it back to the list Cost.

    # # By doing this we can make the input interactive.
    #     Lat=[]
    #     Long=[]
    #     Cost=[]
    #     with open('media/diamond_frame_nodes.csv', newline='') as csvfile:
    #         reader = csv.DictReader(csvfile)
    #         for row in reader:
    #             Lat.append(float(row['Latitude']))
    #             Long.append(float(row['Longitude']))
    #             Cost.append(float(row['Cost']))
        k=0
        # Cost should go as points to be interpolated by bazier curve cost interpolation function.
    
        for i in range(0,N):
            for j in range(0,N):
                Patchz[i,j]=Cost[k]
                k=k+1        
        r=(nTimes-1)*(N-1)+1
        
        zvals=np.zeros((r,r))
        zvals=interpolation(nTimes,N,Patchz) 
        L=float(request.POST['edge_length'])
        L_Diag=2*L*math.sin(alpha)
        Cost_horizontal,Cost_vertical,Cost_diag=Cost_diamomdgraph(N,zvals,nTimes,L,L_Diag)

    # HERE WE ARE TRYING TO PROVIDE A GRAPH FOR THE START NODE TO MITIGATE THE RESTRICTION NEAR SORCE NODE.
    # test the rotation and scaling of patch
    
        nrect=int(request.POST['graph_at_source']) # number of nodes in the graph near source.
    # GRAPH ON RIGHT HAND SIDE
        Patch1x=np.zeros((nrect,nrect))
        Patch1y=np.zeros((nrect,nrect))
        Patch_leftx=np.zeros((nrect,nrect))
        Patch_lefty=np.zeros((nrect,nrect))
        Patch_rightx=np.zeros((nrect,nrect))
        Patch_righty=np.zeros((nrect,nrect))
        

        theta1=math.pi/2+theta-alpha 
        theta2=theta+math.pi+alpha
        Rotation_matrix1=np.array([[math.cos(theta1),math.sin(theta1)],[-math.sin(theta1),math.cos(theta1)]])
        Rotation_matrix2=np.array([[math.cos(theta2),math.sin(theta2)],[-math.sin(theta2),math.cos(theta2)]])

        Patch1x[:,:]=np.array([np.arange(0,nrect,1) for i in range(0,nrect)])
        Patch1y[:,:]=np.array([[j for i in range(0,nrect)] for j in range(0,nrect)])

        l1=float(request.POST['grid_scaling']) 

        Patch1x=Patch1x*l1
        Patch1y=Patch1y*l1
        for i in range(0,nrect):
            for j in range(0,nrect):
                [Patch_leftx[i,j],Patch_lefty[i,j]]=np.dot([Patch1x[i,j],Patch1y[i,j]],Rotation_matrix1) # this is a clocwise rotation
                [Patch_rightx[i,j],Patch_righty[i,j]]=np.dot([Patch1x[i,j],Patch1y[i,j]],Rotation_matrix2)
        # GRAPH ON LEFT HAND SIDE **** 

        Patch_leftx_start=Patch_leftx+Patchx[0,0]
        Patch_lefty_start=Patch_lefty+Patchy[0,0]

        # GRAPH ON RIGHT HAND SIDE  **** 
        
        Patch_rightx_start=Patch_rightx+Patchx[0,0]
        Patch_righty_start=Patch_righty+Patchy[0,0]

     

        Rotation_matrix1=np.array([[math.cos(math.pi/2-alpha),math.sin(math.pi/2-alpha)]\
            ,[-math.sin(math.pi/2-alpha),math.cos(math.pi/2-alpha)]])

        Rotation_matrix2=np.array([[math.cos(alpha),math.sin(alpha)]\
            ,[-math.sin(alpha),math.cos(alpha)]])

        for i in range(0,nrect):
            for j in range(0,nrect):
                [Patch_leftx[i,j],Patch_lefty[i,j]]=np.dot([Patch_leftx[i,j],Patch_lefty[i,j]],Rotation_matrix1) # this is a clocwise rotation
                [Patch_rightx[i,j],Patch_righty[i,j]]=np.dot([Patch_rightx[i,j],Patch_righty[i,j]],Rotation_matrix2)
        

        # GRAPH ON RIGHT HAND SIDE  **** 
        
      
        #  HERE WE ARE PROVIDING THIS SMALL GRAPTH THE INPUTS REQUIRED FOR CALCULATING SHORTEST PATHS NEAR SOURCE.
        #*** LEFTSIDE GRAPH

        Patch_leftx_end=Patch_leftx+Patchx[N-1,N-1]
        Patch_lefty_end=Patch_lefty+Patchy[N-1,N-1]

        Patch_rightx_end=Patch_rightx+Patchx[N-1,N-2]
        Patch_righty_end=Patch_righty+Patchy[N-1,N-2]



        k=0
        Lat_left_start=[1]*nrect**2
        Long_left_start=[1]*nrect**2
        Lat_right_start=[1]*nrect**2
        Long_right_start=[1]*nrect**2

        Lat_left_end=[1]*nrect**2
        Long_left_end=[1]*nrect**2
        Lat_right_end=[1]*nrect**2
        Long_right_end=[1]*nrect**2
        for i in range(0,nrect):
            for j in range(0,nrect):
                Lat_left_start[k]=Patch_lefty_start[i,j]*conversion+start[0]       
                Long_left_start[k]=Patch_leftx_start[i,j]/(111.32*math.cos(Lat_left_start[k]*3.14/180))+start[1]
                Lat_right_start[k]=Patch_righty_start[i,j]*conversion+start[0]       
                Long_right_start[k]=Patch_rightx_start[i,j]/(111.32*math.cos(Lat_right_start[k]*3.14/180))+start[1]

                Lat_left_end[k]=Patch_lefty_end[i,j]*conversion+end[0]       
                Long_left_end[k]=Patch_leftx_end[i,j]/(111.32*math.cos(Lat_left_end[k]*3.14/180))+start[1]
                Lat_right_end[k]=Patch_righty_end[i,j]*conversion+end[0]       
                Long_right_end[k]=Patch_rightx_end[i,j]/(111.32*math.cos(Lat_right_end[k]*3.14/180))+start[1]
                k=k+1
            #return Lat_left_start, Long_left_start,Lat_right_start,Long_right_start,Lat_left_end, Long_left_end,Lat_right_end,Long_right_end
            # with open('E:/Project/Thesis_related_codes/Application_Python/application_shortest_path/shortest_path/media/left_graph.csv', "w") as csv_file:
            #     writer = csv.writer(csv_file, delimiter=',',lineterminator='\n')
            #     writer.writerow(['Latitude_left','Longitude_left','Cost_left','Latitude_right','Longitude_right','Cost_right'])
            #     for i in range(0, len(Lat_left)):
            #         writer.writerow([Lat_left[i],Long_left[i],' ',Lat_right[i],Long_right[i],' '])

        Longitude_left_start=[]
        Latitude_left_start=[]
        Cost_left_start=[]
        Longitude_right_start=[]
        Latitude_right_start=[]
        Cost_right_start=[]

        Longitude_left_end=[]
        Latitude_left_end=[]
        Cost_left_end=[]
        Longitude_right_end=[]
        Latitude_right_end=[]
        Cost_right_end=[]

        if request.method == 'POST':
            Cost_ends = request.FILES['Endgraphs_Cost']
            fs = FileSystemStorage()
            Endgraphs_cost = fs.save(Cost_ends.name, Cost_ends)

        with open('C:/Users/Amit/Desktop/Working/application_shortest_path/shortest_path/media/Endgraphs_cost.csv', newline='') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                Latitude_left_start.append(float(row['Latitude_left_start']))
                Longitude_left_start.append(float(row['Longitude_left_start']))
                Cost_left_start.append(float(row['Cost_left_start']))
                Latitude_right_start.append(float(row['Latitude_right_start']))
                Longitude_right_start.append(float(row['Longitude_right_start']))
                Cost_right_start.append(float(row['Cost_right_start']))

                Latitude_left_end.append(float(row['Latitude_left_end']))
                Longitude_left_end.append(float(row['Longitude_left_end']))
                Cost_left_end.append(float(row['Cost_left_end']))
                Latitude_right_end.append(float(row['Latitude_right_end']))
                Longitude_right_end.append(float(row['Longitude_right_end']))
                Cost_right_end.append(float(row['Cost_right_end']))
        

        # INTERPOLATING COSTS FROM INPUT POINTS
        Cost_left_start=idwr(Latitude_left_start,Longitude_left_start,Cost_left_start,Lat_left_start,Long_left_start)
        Cost_right_start=idwr(Latitude_right_start,Longitude_right_start,Cost_right_start,Lat_right_start,Long_right_start)
        Cost_left_end=idwr(Latitude_left_end,Longitude_left_end,Cost_left_end,Lat_left_end,Long_left_end)
        Cost_right_end=idwr(Latitude_right_end,Longitude_right_end,Cost_right_end,Lat_right_end,Long_right_end)

        # THESE INTERPOLTED COSTS WILL GO INTO THE FINAL INTERPOLATION. HERE THERE WILL BE A
        # DIFFERENT INTERPOLATION FUNCTION IS USED.
        Patch_leftz_start=np.zeros((nrect,nrect))
        Patch_rightz_start=np.zeros((nrect,nrect))
        Patch_rightz_end=np.zeros((nrect,nrect))
        Patch_leftz_end=np.zeros((nrect,nrect))
        k=0 
        for i in range(0,nrect):
            for j in range(0,nrect):
                Patch_leftz_start[i,j]=Cost_left_start[k]
                Patch_rightz_start[i,j]=Cost_right_start[k]
                Patch_leftz_end[i,j]=Cost_left_end[k]
                Patch_rightz_end[i,j]=Cost_right_end[k]
                k=k+1    
            
        r=(nTimes-1)*(nrect-1)+1
        zvals_left_start=np.zeros((r,r))
        zvals_right_start=np.zeros((r,r))
        Costmatrix_left_start=np.zeros((nrect**2,nrect**2)) 
        Costmatrix_right_start=np.zeros((nrect**2,nrect**2))   
        zvals_left_start=interpolation(nTimes,nrect,Patch_leftz_start)
        zvals_right_start=interpolation(nTimes,nrect,Patch_rightz_start)
        # THIS GRAPH IS A COMPLETE GRAPH HENCE THE COST MATRIX IS NOT A SPARSE MATRIX
        Costmatrix_left_start= complete_graph(l1,nrect,zvals_left_start,nTimes)
        Costmatrix_right_start= complete_graph(l1,nrect,zvals_right_start,nTimes)

        zvals_left_end=np.zeros((r,r))
        zvals_right_end=np.zeros((r,r))
        Costmatrix_left_end=np.zeros((nrect**2,nrect**2)) 
        Costmatrix_right_end=np.zeros((nrect**2,nrect**2))   
        zvals_left_end=interpolation(nTimes,nrect,Patch_leftz_end)
        zvals_right_end=interpolation(nTimes,nrect,Patch_rightz_end)
        # THIS GRAPH IS A COMPLETE GRAPH HENCE THE COST MATRIX IS NOT A SPARSE MATRIX
        Costmatrix_left_end= complete_graph(l1,nrect,zvals_left_end,nTimes)
        Costmatrix_right_end= complete_graph(l1,nrect,zvals_right_end,nTimes)

        
        # WE CALCULATE THE SHORTEST PATH NOW USING BELLMAN-FORD ALGORITHM! 

        #%import cost matrix
        start_row=0
        start_col=nrect-1
        end_row=0
        end_col=0
        startnode_start_left=start_row*nrect+start_col
        targetnode_start_left=end_row*nrect+end_col

        for i in range(nrect):            
            Costmatrix_left_start[startnode_start_left,startnode_start_left-i]=1000000 

        sPath_left_start=bellman_ford(nrect**2,startnode_start_left,targetnode_start_left,Costmatrix_left_start,10000000)

        start_row=nrect-1
        start_col=0
        end_row=0
        end_col=0
        startnode_start_right=start_row*nrect+start_col
        targetnode_start_right=end_row*nrect+end_col
        for i in range(nrect):            
            Costmatrix_right_start[startnode_start_right,startnode_start_right-i*nrect]=1000000

        sPath_right_start=bellman_ford(nrect**2,startnode_start_right,targetnode_start_right,Costmatrix_right_start,10000000)

        start_row=nrect-1
        start_col=0
        end_row=0
        end_col=0
        
        startnode_end_left=start_row*nrect+start_col
        targetnode_end_left=end_row*nrect+end_col
        Costmatrix_left_end[startnode_end_left,targetnode_end_left]=1000000
        for i in range(nrect):            
            Costmatrix_left_end[startnode_end_left,startnode_end_left-i*nrect]=1000000 

        sPath_left_end=bellman_ford(nrect**2,startnode_end_left,targetnode_end_left,Costmatrix_left_end,10000000)

        start_row=0
        start_col=0
        end_row=nrect-1
        end_col=0
        startnode_end_right=start_row*nrect+start_col
        targetnode_end_right=end_row*nrect+end_col
        for i in range(nrect):            
            Costmatrix_right_end[startnode_end_right,startnode_end_right-nrect*i]=1000000
        
        sPath_right_end=bellman_ford(nrect**2,startnode_end_right,targetnode_end_right,Costmatrix_right_end,1000000)


        # THIS PART IS FOR THE PURPOSE OF SMALL CORRECTION WE NEED IN THE COST MATRIX.
        # THIS IS TAKING VALUES FROM FRONT END HTML PAGE 

        # A=int(request.POST['edit_row'])
        # B=int(request.POST['edit_column'])
        # D=float(request.POST['cost'])
        # #print(C)
        # Cost_horizontal[A-1,(B-1)]=D 

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
               
        print(sPath)
        #print(Hash)
        # conversion of 2D node arrangement, Patch ---> list of [latitude, longitude].
        
        # for triggering two different submit buttons
        #print(Lat,Long) 
        start= json.dumps(start)
        Lat=json.dumps(Lat)
        Long=json.dumps(Long)
        #Patchx=Patchx.reshape((N)**2,1)
        Patchx=Patchx.tolist()
        Patchx=json.dumps(Patchx)

        #N_K=N_K.reshape(K+1,1)
        N_K=N_K.tolist()
        N_K=json.dumps(N_K)

        #Patchy=Patchy.reshape((N)**2,1)
        Patchy=Patchy.tolist()    
        Patchy=json.dumps(Patchy)

        #Patchz=Patchz.reshape((N)**2,1)
        Patchz=Patchz.tolist()    
        Patchz=json.dumps(Patchz)

        #Patch_rightx_start=Patch_rightx_start.reshape((nrect)**2,1)
        Patch_rightx_start=Patch_rightx_start.tolist()
        Patch_rightx_start=json.dumps(Patch_rightx_start)

        #Patch_righty_start=Patch_righty_start.reshape((nrect)**2,1)
        Patch_righty_start=Patch_righty_start.tolist()
        Patch_righty_start=json.dumps(Patch_righty_start)


        #Patch_leftx_start=Patch_leftx_start.reshape((nrect)**2,1)
        Patch_leftx_start=Patch_leftx_start.tolist()
        Patch_leftx_start=json.dumps(Patch_leftx_start)

        #Patch_lefty_start=Patch_lefty_start.reshape((nrect)**2,1)
        Patch_lefty_start=Patch_lefty_start.tolist()
        Patch_lefty_start=json.dumps(Patch_lefty_start)

        #Patch_rightx_end=Patch_rightx_end.reshape((nrect)**2,1)
        Patch_rightx_end=Patch_rightx_end.tolist()
        Patch_rightx_end=json.dumps(Patch_rightx_end)

        #Patch_righty_end=Patch_righty_end.reshape((nrect)**2,1)
        Patch_righty_end=Patch_righty_end.tolist()
        Patch_righty_end=json.dumps(Patch_righty_end)


        #Patch_leftx_end=Patch_leftx_end.reshape((nrect)**2,1)
        Patch_leftx_end=Patch_leftx_end.tolist()
        Patch_leftx_end=json.dumps(Patch_leftx_end)

        #Patch_lefty_end=Patch_lefty_end.reshape((nrect)**2,1)
        Patch_lefty_end=Patch_lefty_end.tolist()
        Patch_lefty_end=json.dumps(Patch_lefty_end)

        #sPath=sPath.reshape((sPath.shape[0]*3),1)
        sPath=sPath.tolist()
        sPath=json.dumps(sPath)

        sPath_left_start=json.dumps(sPath_left_start)
        sPath_right_start=json.dumps(sPath_right_start)

        #Hash=Hash.reshape((n*3),1)
        Hash=Hash.tolist()
        Hash=json.dumps(Hash)
        
        #Cost_horizontal=Cost_horizontal.reshape((N)**2,1)
        Cost_horizontal=Cost_horizontal.tolist()
        Cost_horizontal=json.dumps(Cost_horizontal)
        #Cost_vertical=Cost_vertical.reshape((N)**2,1)
        Cost_vertical=Cost_vertical.tolist()
        Cost_vertical=json.dumps(Cost_vertical)
        #Cost_diag=Cost_diag.reshape((N)**2,1)
        Cost_diag=Cost_diag.tolist()
        Cost_diag=json.dumps(Cost_diag)
        N=json.dumps(N)
        context={'N':N}
        nrect=json.dumps(nrect)
        context={'nrect':nrect}
        K=json.dumps(K)
        # context={'K':K}
        # context = {'sPath': sPath}
        # context = {'sPath_left': sPath_left}
        # context = {'sPath_right': sPath_right}
        # context = {'Hash': Hash}
        # context = {'Cost_horizontal': Cost_horizontal}
        # context = {'Patchx': Patchx}
        # context = {'Patchy': Patchy}
        # context = {'Patch_rightx': Patch_rightx}
        # context = {'Patch_righty': Patch_righty}
        # context = {'Patch_leftx': Patch_leftx}
        # context = {'Patch_lefty': Patch_lefty}
        # context = {'Cost_vertical': Cost_vertical}
        # context = {'Cost_diag': Cost_diag}
        # context = {'start': start}
        # context = {'Long': Long}
        # context = {'Lat': Lat}
        request.session['N_K']=N_K
        request.session['nrect']=nrect
        request.session['N'] = N
        request.session['K'] = K
        request.session['Hash'] = Hash
        request.session['Cost_horizontal'] = Cost_horizontal
        request.session['Cost_vertical'] = Cost_vertical
        request.session['Cost_diag'] = Cost_diag
        request.session['nrect'] = nrect
        request.session['start'] = start
        request.session['Patchx'] = Patchx
        request.session['Patchy'] = Patchy
        request.session['Lat'] = Lat
        request.session['Long'] = Long
        request.session['Patch_leftx_start'] = Patch_leftx_start
        request.session['Patch_rightx_start'] = Patch_rightx_start
        request.session['Patch_leftx_end'] = Patch_leftx_end
        request.session['Patch_rightx_end'] = Patch_rightx_end
        request.session['Patch_lefty_start'] = Patch_lefty_start
        request.session['Patch_righty_start'] = Patch_righty_start
        request.session['Patch_lefty_end'] = Patch_lefty_end
        request.session['Patch_righty_end'] = Patch_righty_end
        request.session['sPath_left_start'] = sPath_left_start
        request.session['sPath_right_start'] = sPath_right_start
        request.session['sPath_left_end'] = sPath_left_end
        request.session['sPath_right_end'] = sPath_right_end
        request.session['sPath'] = sPath
    #print(Hash)
    
    # request.session['Patchz'] = Patchz
    #print(1)
    return render(request,"shortest_path1.html",{'sPath':sPath,'sPath_left_start':sPath_left_start,'sPath_right_start':sPath_right_start,\
        'sPath_left_end':sPath_left_end,'sPath_right_end':sPath_right_end,\
    'N':N,'K':K,'Cost_horizontal':Cost_horizontal,'Cost_vertical':Cost_vertical,'Hash':Hash,'Cost_diag':Cost_diag,\
    'nrect':nrect,'start':start,'Patchx':Patchx,'Lat':Lat,'Long':Long,'Patchy':Patchy,\
    'Patch_rightx_start':Patch_rightx_start,'Patch_leftx_start':Patch_leftx_start,\
         'Patch_righty_start':Patch_righty_start,'Patch_lefty_start':Patch_lefty_start, 'Patch_rightx_end':Patch_rightx_end,'Patch_leftx_end':Patch_leftx_end,\
         'Patch_righty_end':Patch_righty_end,'Patch_lefty_end':Patch_lefty_end})



def Modify_cost(request):
    sPath=request.session['sPath']
    
    Patchx=request.session['Patchx']

    Patchy=request.session['Patchy']
  
    Lat=request.session['Lat']
    Long=request.session['Long'] 
    Patch_leftx_start=request.session['Patch_leftx_start'] 
    Patch_rightx_start=request.session['Patch_rightx_start']
    Patch_lefty_start=request.session['Patch_lefty_start']
    Patch_righty_start=request.session['Patch_righty_start']
    sPath_left_start=request.session['sPath_left_start']
    sPath_right_start=request.session['sPath_right_start']

    Patch_leftx_end=request.session['Patch_leftx_end'] 
    Patch_rightx_end=request.session['Patch_rightx_end']
    Patch_lefty_end=request.session['Patch_lefty_end']
    Patch_righty_end=request.session['Patch_righty_end']
    sPath_left_end=request.session['sPath_left_end']
    sPath_right_end=request.session['sPath_right_end']
    
    nrect=request.session['nrect']
    conversion=1/111
    N=request.session['N']
    N=json.loads(N)
    N_K=request.session['N_K']
    N_K=json.loads(N_K)
    N_K=np.array(N_K)
    K=request.session['K']
    K=json.loads(K)
    n=int(1/48*((K+2)*(K+4)*(2*K+6)))
    Hash=request.session['Hash'] 
    
    Cost_horizontal=request.session['Cost_horizontal']    
    Cost_vertical= request.session['Cost_vertical']
    Cost_diag= request.session['Cost_diag']
    Cost_horizontal=json.loads(Cost_horizontal)
    Cost_vertical=json.loads(Cost_vertical)
    Cost_diag=json.loads(Cost_diag)
    Cost_horizontal=np.array(Cost_horizontal)
    Cost_vertical=np.array(Cost_vertical)
    Cost_diag=np.array(Cost_diag)
    Hash=json.loads(Hash)
    Hash=np.array(Hash) 
    #nrect= request.session['nrect'] 
    start= request.session['start'] 
    #Patchx= request.session['Patchx'] 
    #Patchy= request.session['Patchy']
    
    
    iD=int(request.POST['id_nodes'])
    #print(type(iD))
    Cost_horizontal[int(iD//N),int(iD%N)+1]=float(request.POST['Cost_horizontal'])
    #print(Cost_horizontal[int(iD//N),int(iD%N)])
    Cost_vertical[int(iD//N)+1,int(iD%N)]=float(request.POST['Cost_vertical'])
    Cost_diag[int(iD//N)+1,int(iD%N)+1]=float(request.POST['Cost_diag'])


    #print(Cost_diag[:,:])
    #print(end_time-start_time)
    sPath1=np.zeros((K,3))
    Node_val=node_val(Hash,Cost_horizontal,Cost_vertical,Cost_diag,N_K,K)
    sPath1[K-1,:]=Node_val[n-1,:]
    
    for m in range(K-2,-1,-1):
        sPath1[m,:]=Node_val[int(sPath1[m+1,1]),:]

    sPath1=sPath1.tolist()
    sPath1=json.dumps(sPath1)
    
    context = {'sPath1': sPath1}
    Cost_vertical=Cost_vertical.tolist()
    Cost_vertical=json.dumps(Cost_vertical)
    Cost_horizontal=Cost_horizontal.tolist()
    Cost_horizontal=json.dumps(Cost_horizontal)
    Cost_diag=Cost_diag.tolist()
    Cost_diag=json.dumps(Cost_diag)
    N_K=N_K.tolist()
    N_K=json.dumps(N_K)
    N=json.dumps(N)
   

    #Patchy=Patchy.reshape((N)**2,1)


 

    #Hash=Hash.reshape((n*3),1)
    Hash=Hash.tolist()
    Hash=json.dumps(Hash)

    # start= json.dumps(start)
    # Lat=json.dumps(Lat)
    # Long=json.dumps(Long)


    
    return render(request,"Modify_cost.html",{'sPath1':sPath1,'sPath':sPath,'sPath_left_start':sPath_left_start,'sPath_right_start':sPath_right_start,\
        'sPath_left_end':sPath_left_end,'sPath_right_end':sPath_right_end,\
            'N':N,'K':K,'Cost_horizontal':Cost_horizontal,'Cost_vertical':Cost_vertical,'Hash':Hash,'Cost_diag':Cost_diag,
    'nrect':nrect,'start':start,'Patchx':Patchx,'Lat':Lat,'Long':Long,'Patchy':Patchy,
    'Patch_rightx_start':Patch_rightx_start,'Patch_leftx_start':Patch_leftx_start,\
         'Patch_righty_start':Patch_righty_start,'Patch_lefty_start':Patch_lefty_start, 'Patch_rightx_end':Patch_rightx_end,'Patch_leftx_end':Patch_leftx_end,\
         'Patch_righty_end':Patch_righty_end,'Patch_lefty_end':Patch_lefty_end})

    #return redirect("shortest_path1.html",{'sPath1':sPath1})

    #print(sPath)

#     gmap5 = gmplot.GoogleMapPlotter(start[0],start[1],13) #assign a map from Gmapplotter 
#     gmap5.apikey = request.POST['enter your gmap API key']    
#     colors = ['red','blue','green','purple','orange','yellow','pink','white']

    # for i in range(0,N-1):
    #     for j in range(0,N-1):           
    #         lat1=Patchy[i,j]*conversion+start[0]   
    #         long1=Patchx[i,j]/(111.32*math.cos(lat1*3.14/180))+start[1]
    #         lat2=Patchy[i,j+1]*conversion+start[0]
    #         long2=Patchx[i,j+1]/(111.32*math.cos(lat2*3.14/180))+start[1]
    #         gmap5.plot([lat1,lat2],[long1,long2], colors[2], edge_width = 1.5)
                    
    #         lat2=Patchy[i+1,j]*conversion+start[0] 
    #         long2=Patchx[i+1,j]/(111.32*math.cos(lat2*3.14/180))+start[1]
                
    #         gmap5.plot([lat1,lat2], [long1,long2], colors[3], edge_width = 1.5)

    # ends=N-1
    # for i in range(0,N-1):
    #     lat1=Patchy[i,ends]*conversion+start[0]   
    #     long1=Patchx[i,ends]/(111.32*math.cos(lat1*3.14/180))+start[1]
    #     lat2=Patchy[i+1,ends]*conversion+start[0]   
    #     long2=Patchx[i+1,ends]/(111.32*math.cos(lat2*3.14/180))+start[1]
    #     gmap5.plot([lat1,lat2], [long1,long2], colors[2], edge_width = 1.5)

    #     lat1=Patchy[ends,i]*conversion+start[0]   
    #     long1=Patchx[ends,i]/(111.32*math.cos(lat1*3.14/180))+start[1]
    #     lat2=Patchy[ends,i+1]*conversion+start[0]   
    #     long2=Patchx[ends,i+1]/(111.32*math.cos(lat2*3.14/180))+start[1]
    #     gmap5.plot([lat1,lat2], [long1,long2], colors[3], edge_width = 1.5)    
    
    # for i in range(0,nrect-1):
    #     for j in range(0,nrect-1):           
    #         lat1=Patch_righty[i,j]*conversion+start[0]   
    #         long1=Patch_rightx[i,j]/(111.32*math.cos(lat1*3.14/180))+start[1]
    #         lat2=Patch_righty[i,j+1]*conversion+start[0]
    #         long2=Patch_rightx[i,j+1]/(111.32*math.cos(lat2*3.14/180))+start[1]
    #         gmap5.plot([lat1,lat2],[long1,long2], colors[2], edge_width = 1.5)
                    
    #         lat2=Patch_righty[i+1,j]*conversion+start[0] 
    #         long2=Patch_rightx[i+1,j]/(111.32*math.cos(lat2*3.14/180))+start[1]                   
    #         gmap5.plot([lat1,lat2], [long1,long2], colors[3], edge_width = 1.5)

    # ends=nrect-1
    # for i in range(0,nrect-1):
    #     lat1=Patch_righty[i,ends]*conversion+start[0]   
    #     long1=Patch_rightx[i,ends]/(111.32*math.cos(lat1*3.14/180))+start[1]
    #     lat2=Patch_righty[i+1,ends]*conversion+start[0]   
    #     long2=Patch_rightx[i+1,ends]/(111.32*math.cos(lat2*3.14/180))+start[1]
    #     gmap5.plot([lat1,lat2], [long1,long2], colors[2], edge_width = 1.5)

    #     lat1=Patch_righty[ends,i]*conversion+start[0]   
    #     long1=Patch_rightx[ends,i]/(111.32*math.cos(lat1*3.14/180))+start[1]
    #     lat2=Patch_righty[ends,i+1]*conversion+start[0]   
    #     long2=Patch_rightx[ends,i+1]/(111.32*math.cos(lat2*3.14/180))+start[1]
    #     gmap5.plot([lat1,lat2], [long1,long2], colors[3], edge_width = 1.5)

    # for i in range(0,nrect-1):
    #     for j in range(0,nrect-1):           
    #         lat1=Patch_lefty[i,j]*conversion+start[0]   
    #         long1=Patch_leftx[i,j]/(111.32*math.cos(lat1*3.14/180))+start[1]
    #         lat2=Patch_lefty[i,j+1]*conversion+start[0]
    #         long2=Patch_leftx[i,j+1]/(111.32*math.cos(lat2*3.14/180))+start[1]
    #         gmap5.plot([lat1,lat2],[long1,long2], colors[2], edge_width = 1.5)
                    
    #         lat2=Patch_lefty[i+1,j]*conversion+start[0] 
    #         long2=Patch_leftx[i+1,j]/(111.32*math.cos(lat2*3.14/180))+start[1]                   
    #         gmap5.plot([lat1,lat2], [long1,long2], colors[3], edge_width = 1.5)

    # ends=nrect-1
    # for i in range(0,nrect-1):
    #     lat1=Patch_lefty[i,ends]*conversion+start[0]   
    #     long1=Patch_leftx[i,ends]/(111.32*math.cos(lat1*3.14/180))+start[1]
    #     lat2=Patch_lefty[i+1,ends]*conversion+start[0]   
    #     long2=Patch_leftx[i+1,ends]/(111.32*math.cos(lat2*3.14/180))+start[1]
    #     gmap5.plot([lat1,lat2], [long1,long2], colors[2], edge_width = 1.5)

    #     lat1=Patch_lefty[ends,i]*conversion+start[0]   
    #     long1=Patch_leftx[ends,i]/(111.32*math.cos(lat1*3.14/180))+start[1]
    #     lat2=Patch_lefty[ends,i+1]*conversion+start[0]   
    #     long2=Patch_leftx[ends,i+1]/(111.32*math.cos(lat2*3.14/180))+start[1]
    #     gmap5.plot([lat1,lat2], [long1,long2], colors[3], edge_width = 1.5)
        
    # #plotting the shortest route on 2D                    
    # gmap5.draw( "Template/map.html" )
    # return render(request,"map.html")
#    for v1 in range(0,K):        
      
#         if (sPath[v1,0]>0.0):
#             [N1,N2] = [int(sPath[v1,1]),int(sPath[v1,2])]
#             [j1,i11,i12]=[Hash[N1,0],Hash[N1,1],Hash[N1,2]]
#             [j2,i21,i22]=[Hash[N2,0],Hash[N2,1],Hash[N2,2]]
#             [X1,Y1,X2,Y2,X3,Y3,X4,Y4]=[int((j1-i11)/2),int((j1+i11)/2),int((j1-i12)/2),int((j1+i12)/2),int((j2-i21)/2),int((j2+i21)/2),int((j2-i22)/2),int((j2+i22)/2)]
            
#             lat1 = Patchy[X1,Y1]*conversion+start[0]
#             lat2 = Patchy[X3,Y3]*conversion+start[0]            
#             long1 = Patchx[X1,Y1]/(111.32*math.cos(lat1*3.14/180))+start[1]
#             long2 = Patchx[X3,Y3]/(111.32*math.cos(lat2*3.14/180))+start[1]
#             gmap5.plot([lat1, lat2],[long1, long2],colors[2], edge_width = 2)

#             lat1 = Patchy[X2,Y2]*conversion+start[0]
#             lat2 = Patchy[X4,Y4]*conversion+start[0]            
#             long1 = Patchx[X2,Y2]/(111.32*math.cos(lat1*3.14/180))+start[1]
#             long2 = Patchx[X4,Y4]/(111.32*math.cos(lat2*3.14/180))+start[1]
#             gmap5.plot([lat1, lat2],[long1, long2],colors[3], edge_width = 2)
#     k=0
#     for v1 in range(0,len(sPath_left)-1):   
#         [N1,N2] = [int(sPath_left[v1]),int(sPath_left[v1+1])]
#         X1=N1//nrect
#         X2=N2//nrect
#         Y1=N1%nrect
#         Y2=N2%nrect       
#         lat1 = Patch_lefty[X1,Y1]*conversion+start[0]
#         lat2 = Patch_lefty[X2,Y2]*conversion+start[0]            
#         long1 = Patch_leftx[X1,Y1]/(111.32*math.cos(lat1*3.14/180))+start[1]
#         long2 = Patch_leftx[X2,Y2]/(111.32*math.cos(lat2*3.14/180))+start[1]
#         gmap5.plot([lat1, lat2],[long1, long2],colors[2], edge_width = 2)
#     k=0
#     for v1 in range(0,len(sPath_right)-1):
#         [N1,N2] = [int(sPath_right[v1]),int(sPath_right[v1+1])]
#         X1=N1//nrect
#         X2=N2//nrect
#         Y1=N1%nrect
#         Y2=N2%nrect       
#         lat1 = Patch_righty[X1,Y1]*conversion+start[0]
#         lat2 = Patch_righty[X2,Y2]*conversion+start[0]            
#         long1 = Patch_rightx[X1,Y1]/(111.32*math.cos(lat1*3.14/180))+start[1]
#         long2 = Patch_rightx[X2,Y2]/(111.32*math.cos(lat2*3.14/180))+start[1]
    #     gmap5.plot([lat1, lat2],[long1, long2],colors[3], edge_width = 2)


