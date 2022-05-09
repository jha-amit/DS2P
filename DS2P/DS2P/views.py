from django.shortcuts import redirect, render
import sys
import numpy as np
import random
import math
import gmplot
import csv
from cython_node_val import node_val
from Cost_modified import interpolation,Cost_diamondgraph,Patchpoints,hash1,complete_graph,Patch_terminal,interpolation_radial,Cost_radial
import pickle
from .forms import UploadFileForm
from django.conf import settings
from django.core.files.storage import FileSystemStorage
from django.core.files.storage import FileSystemStorage
from django.db import models
from numpy import genfromtxt
import json 
from django.core.files.storage import default_storage
from IDW import IDW
from bellman_ford import bellman_ford
from django.views.decorators.csrf import csrf_exempt
from django.http import HttpResponse, JsonResponse
from django.core import serializers
import time

val=None


def home(request):
   return render(request,'Home_test.html',{"name":"user"}) 

def plot(request):
    t1=time.time()
    if  request.POST:
        start_lat=float(request.POST['Start_point_latitude']) 
        start_long=float(request.POST['Start_point_longitude'])
        end_lat=float(request.POST['End_point_latitude'])  
        end_long=float(request.POST['End_point_longitude'])
        d=int(request.POST['distance'])
        K=int(request.POST['number_of_layers'])
        alpha=float(request.POST['edge_slope'])
        gmap_license=request.POST["enter_your_gmap_API_key"]

    conversion=1/111 #conversion from change in latitude to km.
                    #conversion for km to latitude = (111.32*math.cos(lat[k]*math.pi/180)) 
    
    end=[end_lat,end_long]   
    start=[start_lat,start_long]
    gmap_license="https://maps.googleapis.com/maps/api/js?key="+gmap_license
    #print(gmap_license)
    if ((K%2!=0) or (K < 0)):
        return render(request,"Home.html",{'shortest_path1':"Enter Correct Inputs and Press Input button"})
    else:
        # This parameter is to control the actual edge length. it divides 1 unit in l parts so we
                #can shrink or expand the grid for our uses. here i am using it for optimal size of grid to acheiv
                #  fairly accurate cost interpolations. 
        N=int(K/2)+1
        Patchx=np.zeros((N,N))
        Patchy=np.zeros((N,N))
        Patchz=np.zeros((N,N))
        
        # 111 km per degree latitude change and 111.17km for per degree longitude change 
        #theta=-math.pi/3

        theta=math.atan((((end_lat-start_lat)*111)/((end_long-start_long)*111.2))) # need correction for eliptical co-ordinates
        
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
           
            
        d_half=int(d/2)
        
        
        def distance(x1,y1,x2,y2):
            dist=math.sqrt((x1-x2)**2+(y1-y2)**2)
            return dist

        delta=distance(Patchx[0,1],Patchy[0,1],Patchx[1,0],Patchy[1,0])
    

        incremental_radius=float(request.POST['increment'])
        Center_to_diamondtip=float(request.POST["dist_center_to_tip"])
        
            
        n_circum_start=int(request.POST['n_circum_start'])
        n_circum_end=int(request.POST['n_circum_end'])

       

        centery_start= Patchy[0,0]+ Center_to_diamondtip*math.sin(theta+(math.pi))

        centerx_start= Patchx[0,0]+ Center_to_diamondtip*math.cos(theta+(math.pi))
        
        center_start_lat = start[0] + centery_start*conversion

        center_start_long= start[1] +centerx_start/(111.32*math.cos(center_start_lat*3.14/180))

        centery_end= Patchy[N-1,N-1]+ Center_to_diamondtip*math.sin(theta)

        centerx_end= Patchx[N-1,N-1]+ Center_to_diamondtip*math.cos(theta)
        
        center_end_lat = start[0] + centery_end*conversion

        center_end_long= start[1] +centerx_end/(111.32*math.cos(center_end_lat*3.14/180))

        user_defined_radius_start = distance(Patchx[d_half,0],Patchy[d_half,0],centerx_start,centery_start)
        user_defined_radius_end = distance(Patchx[N-1-d_half,N-1],Patchy[N-1-d_half,N-1],centerx_end,centery_end)
        #print(user_defined_radius_start)

        n_radial_end=int(user_defined_radius_end/incremental_radius)+1
        n_radial_start=int(user_defined_radius_start/incremental_radius)+1
        

        radial_nodes_lat_start, radial_nodes_long_start,\
            radial_nodes_starty, radial_nodes_startx,radial_nodes_lat_end, radial_nodes_long_end,\
            radial_nodes_endy, radial_nodes_endx\
                     = Patch_terminal(center_start_lat,center_start_long,center_end_lat,center_end_long,\
                          n_circum_start, n_radial_start,theta,user_defined_radius_start, n_circum_end,\
                               n_radial_end,user_defined_radius_end,conversion)

        # radial_nodes_lat_end, radial_nodes_long_end,\
        #     radial_nodes_endy, radial_nodes_endx\
        #              = Patch_terminal(center_end_lat,\
        #                     center_end_long,n_circum_end,n_radial_end,theta,\
        #                         user_defined_radius_end,conversion)

        radial_nodes_lat_start=radial_nodes_lat_start.tolist()
        #radial_nodes_lat_start=json.dumps(radial_nodes_lat_start)

        radial_nodes_long_start=radial_nodes_long_start.tolist()
        #radial_nodes_long_start=json.dumps(radial_nodes_long_start)

        radial_nodes_starty=radial_nodes_starty.tolist()
        #radial_nodes_starty=json.dumps(radial_nodes_starty)

        radial_nodes_startx=radial_nodes_startx.tolist()
        #radial_nodes_startx=json.dumps(radial_nodes_startx)

        radial_nodes_lat_end=radial_nodes_lat_end.tolist()
        #radial_nodes_lat_end=json.dumps(radial_nodes_lat_end)

        radial_nodes_long_end=radial_nodes_long_end.tolist()
        #radial_nodes_long_end=json.dumps(radial_nodes_long_end)

        radial_nodes_endy=radial_nodes_endy.tolist()
        #radial_nodes_endy=json.dumps(radial_nodes_endy)

        radial_nodes_endx=radial_nodes_endx.tolist()
        #radial_nodes_endx=json.dumps(radial_nodes_endx)

        Patchx=Patchx.tolist()
        #Patchx=json.dumps(Patchx)
        Patchy=Patchy.tolist()    
        #Patchy=json.dumps(Patchy)    
            
        #Patchz=json.dumps(Patchz)

        # n_circum_start=json.dumps(n_circum_start)
        # n_circum_end=json.dumps(n_circum_end)
        # n_radial_start=json.dumps(n_radial_start)
        # n_radial_end=json.dumps(n_radial_end)

        #gmap_license=json.dumps(gmap_license)
        #n_radial_start=json.dumps(n_radial_start)
        request.session['n_radial_start']=n_radial_start
        request.session['n_radial_end']=n_radial_end
        request.session['n_circum_start']=n_circum_start
        request.session['n_circum_end']=n_circum_end

        #start=json.dumps(start) 
        # d=json.dumps(d)
        # K=json.dumps(K)
        # N=json.dumps(N)
        # Lat=json.dumps(Lat)
        # Long=json.dumps(Long)
        request.session['d'] = d

        request.session['N'] = N
        request.session['alpha'] = alpha
        request.session['gmap_license']=gmap_license
        request.session['K'] = K     
        request.session['start'] = start
        request.session['Patchx'] = Patchx
        request.session['Patchy'] = Patchy
        request.session['Lat'] = Lat
        request.session['Long'] = Long
        request.session['user_defined_radius_start']=user_defined_radius_start
        request.session['user_defined_radius_end']=user_defined_radius_end

        request.session['radial_nodes_long_start']=radial_nodes_long_start
        request.session['radial_nodes_lat_start']=radial_nodes_lat_start
       
        request.session['radial_nodes_starty']=radial_nodes_starty
        request.session['radial_nodes_startx']=radial_nodes_startx
        request.session['radial_nodes_lat_end']=radial_nodes_lat_end
        request.session['radial_nodes_long_end']=radial_nodes_long_end
        request.session['radial_nodes_endy']=radial_nodes_endy
        request.session['radial_nodes_endx']=radial_nodes_endx  

        t2=time.time()
        print(t2-t1)
    return JsonResponse ({'radial_nodes_long_start' : radial_nodes_long_start,\
            'radial_nodes_lat_start': radial_nodes_lat_start,'radial_nodes_long_end' : radial_nodes_long_end,\
            'radial_nodes_lat_end': radial_nodes_lat_end,'N':N,'K':K,'d':d,'start':start,'Patchx':Patchx,'Lat':Lat,\
        'Long':Long,'Patchy':Patchy,'gmap_license':gmap_license,'n_circum_start':n_circum_start,\
            'n_radial_start':n_radial_start,'n_circum_end':n_circum_end,'n_radial_end':n_radial_end})

def cost_manipulate(request):
    Lat=request.session['Lat']
    Long=request.session['Long']
    alpha=request.session['alpha']
    N=request.session['N']
    d=request.session['d']
    K=request.session['K']
    Lat=request.session['Lat']
    Long=request.session['Long']
    Patchx=request.session['Patchx']
    Patchy=request.session['Patchy']   
    start=request.session['start']
    Patchz=np.zeros((N,N))
    
    Patchx=np.array(Patchx)
    Patchy=np.array(Patchy)
# surface patch density and edgelength
    if  request.is_ajax() and request.POST:
        nTimes=int(request.POST['Grid_density'])
        L=float(request.POST['edge_length'])

        
    
    L_Diag=2*L*math.sin(alpha) 
    N_j= N # in case the graph is not square
   
    Longitude_input=[]
    Latitude_input=[]
    Cost_input=[]
    
    myfile = request.FILES['Cost_matrix']
    fs = FileSystemStorage()
    Cost_matrix = fs.save(myfile.name, myfile)

    # Reading the file from filestorage 
    with open('media/Cost_matrix.csv', newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            Latitude_input.append(float(row['Latitude']))
            Longitude_input.append(float(row['Longitude']))
            Cost_input.append(float(row['Cost']))
                
                        
    # Inverse distance weight algorithm. We get a list of interpolated points for the unknown xi, yi,
    # using known points x, y, z. 
    Latitude_input=np.array(Latitude_input)
    Longitude_input=np.array(Longitude_input)
    Cost_input=np.array(Cost_input)
    Lat=np.array(Lat)
    Long=np.array(Long)

    #Interpolate costs at nodes

    Cost=IDW(Latitude_input,Longitude_input,Cost_input,Lat,Long,1,1)

    # receive list of nodes whose value is infinity
    node_ids= request.POST.getlist('nodeid_forcostchange[]')

    print(node_ids)

    for i in node_ids:
        Cost[int(i)]=10000000000    
        print(int(i))
    # Write the diamond node costs to a .csv file to enable user for editing  
    with open('media/diamond_frame_nodes.csv', "w") as csv_file:
        writer = csv.writer(csv_file, delimiter=',',lineterminator='\n')
        writer.writerow(['Latitude','Longitude','Cost'])
        for i in range(0, len(Lat)):
            writer.writerow([Lat[i],Long[i],Cost[i]])

    k=0    
    for i in range(0,N):
        for j in range(0,N):
            Patchz[i,j]=Cost[k]
            k=k+1        
    r=(nTimes-1)*(N-1)+1
    
    zvals=np.zeros((r,r))
    zvals=interpolation(nTimes,N,Patchz)

    
    print(L)
    Cost_horizontal,Cost_vertical,Cost_diag=Cost_diamondgraph(N,zvals,nTimes,L,L_Diag)
    print(Cost_horizontal[4,8])
    with open('media/cost_horizontal.csv', "w") as csv_file:
        writer = csv.writer(csv_file, delimiter=',',lineterminator='\n')
        #writer.writerow(['Latitude','Longitude','Cost'])
        for i in range(0, N):
            writer.writerow(Cost_horizontal[i,:])
    #Cost_horizontal=np.genfromtxt("Cost_horizontal_right _8.csv",delimiter=",")
    #Cost_vertical=np.genfromtxt("Cost_vertical_down _8.csv", delimiter=",")    
    Cost_horizontal=Cost_horizontal.tolist()
    Cost_vertical=Cost_vertical.tolist()
    Cost_diag=Cost_diag.tolist()
    
    request.session['K']=K
    request.session['d']=d
    request.session['L']=L
    
    # request.session['Cost_horizontal'] = Cost_horizontal
    # request.session['Cost_vertical'] = Cost_vertical
    # request.session['Cost_diag'] = Cost_diag
    request.session['Cost_horizontal'] = Cost_horizontal
    request.session['Cost_vertical'] = Cost_vertical
    request.session['Cost_diag'] = Cost_diag
    request.session['N']=N

    return JsonResponse({'N':N})

def radial_SP(request):
    t1=time.time()
      
    n_radial_start=request.session['n_radial_start']
    n_radial_end=request.session['n_radial_end']
    n_circum_start=request.session['n_circum_start']
    n_circum_end=request.session['n_circum_end']

    radial_nodes_lat_start=request.session['radial_nodes_lat_start']
    radial_nodes_lat_start=np.array(radial_nodes_lat_start)

    radial_nodes_long_start=request.session['radial_nodes_long_start']
    radial_nodes_long_start=np.array(radial_nodes_long_start)

    radial_nodes_starty=request.session['radial_nodes_starty']
    radial_nodes_starty=np.array(radial_nodes_starty)

    radial_nodes_startx=request.session['radial_nodes_startx']
    radial_nodes_startx=np.array(radial_nodes_startx)

    
    radial_nodes_lat_end=request.session['radial_nodes_lat_end']
    radial_nodes_lat_end=np.array(radial_nodes_lat_end)

    radial_nodes_long_end=request.session['radial_nodes_long_end']
    radial_nodes_long_end=np.array(radial_nodes_long_end)

    radial_nodes_endy=request.session['radial_nodes_endy']
    radial_nodes_endy=np.array(radial_nodes_endy)

    radial_nodes_endx=request.session['radial_nodes_endx']
    radial_nodes_endx=np.array(radial_nodes_endx)

    alpha=request.session['alpha']
    
    user_defined_radius_start=request.session['user_defined_radius_start']
    user_defined_radius_end=request.session['user_defined_radius_end']

    d=request.session['d']
    K=request.session['K']
    N=request.session['N']
    start=request.session['start']

    #start=json.loads(start)
    Longitude_start=[]
    Latitude_start=[]
    Cost_start=[]

    Longitude_end=[]
    Latitude_end=[]
    Cost_end=[]

    #if  request.is_ajax() and request.POST:
    nTimes=int(request.POST['Grid_density'])

    myfile = request.FILES["Endgraphs_Cost_start"]
    fs = FileSystemStorage()
    Endgraphs_Cost_start = fs.save(myfile.name, myfile)

    myfile = request.FILES["Endgraphs_Cost_end"]
    fs = FileSystemStorage()
    Endgraphs_Cost_end = fs.save(myfile.name, myfile)


    # Reading the file from filestorage 
    with open('media/radial_start.csv', newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            Latitude_start.append(float(row['Latitude']))
            Longitude_start.append(float(row['Longitude']))
            Cost_start.append(float(row['Cost']))

    with open('media/radial_end.csv', newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            Latitude_end.append(float(row['Latitude']))
            Longitude_end.append(float(row['Longitude']))
            Cost_end.append(float(row['Cost']))
            
    Cost_start_input=IDW(np.array(Latitude_start),np.array(Longitude_start),np.array(Cost_start),radial_nodes_lat_start.\
        reshape(int(n_circum_start)*n_radial_start,), radial_nodes_long_start.reshape((n_circum_start)*n_radial_start,),1,1)

    

    Cost_end_input=IDW(np.array(Latitude_end),np.array(Longitude_end),np.array(Cost_end),radial_nodes_lat_end.\
        reshape(int(n_circum_end)*n_radial_end,), radial_nodes_long_end.reshape((n_circum_end)*n_radial_end,),1,1)

    
    Patchz_start=np.array(Cost_start_input).reshape(n_radial_start,n_circum_start)

    Patchz_end=np.array(Cost_end_input).reshape(n_radial_end,n_circum_end)            

    r_start_circum=(nTimes-1)*(n_circum_start-1)+1
    r_start_radial=(nTimes-1)*(n_radial_start-1)+1

    r_end_circum=(nTimes-1)*(n_circum_end-1)+1
    r_end_radial=(nTimes-1)*(n_radial_end-1)+1

    # r_end_circum=(nTimes-1)*(n_circum_end)+1
    # r_end_radial=(nTimes-1)*(n_radial_start-1)+1

    zvals_start=np.zeros((r_start_radial,r_start_circum))
    zvals_end=np.zeros((r_end_radial,r_end_circum))

    
    xvals_start=np.zeros((r_start_radial,r_start_circum))
    yvals_start=np.zeros((r_start_radial,r_start_circum))

    

    zvals_start= interpolation_radial(nTimes,zvals_start,Patchz_start)
    zvals_end= interpolation_radial(nTimes,zvals_end,Patchz_end) 
    


    Cost_radial_start= Cost_radial(n_circum_start,zvals_start,nTimes,n_radial_start,\
                                radial_nodes_startx, radial_nodes_starty)

    #print(radial_nodes_endx.shape,radial_nodes_endx,n_circum_end,n_radial_end)

    Cost_radial_end= Cost_radial(n_circum_end,zvals_end,nTimes,n_radial_end,\
                                radial_nodes_endx, radial_nodes_endy)


    N_start=int(n_circum_start*n_radial_start)
    N_end=int(n_circum_end*n_radial_end)
    
    start_node_start=int(1.0)
    start_node_end=int(1.0)
    # get the node corresponding to the angle subtained by secant of width 'd'. 
    first_target_node_start = math.ceil(math.asin(d/(2* user_defined_radius_start))*(n_circum_start-1)/(2*math.pi))
    first_target_node_end = math.ceil(math.asin(d/(2* user_defined_radius_end))*(n_circum_end-1)/(2*math.pi))


    target_nodes_start1= [(n_circum_start*(n_radial_start-1)+int(i)) for i in np.linspace(0,first_target_node_start,first_target_node_start+1)]
    target_nodes_start2= [(n_circum_start*(n_radial_start)-1-int(i)) for i in np.linspace(0,first_target_node_start,first_target_node_start+1)]
    target_nodes_start= target_nodes_start1+target_nodes_start2
    
    #print(target_nodes_start)
    
    target_nodes_end1= [(n_circum_end*(n_radial_end-1)+int(i)) for i in np.linspace(0,first_target_node_end,first_target_node_end+1)]
    target_nodes_end2= [(n_circum_end*(n_radial_end)-1-int(i)) for i in np.linspace(0,first_target_node_end,first_target_node_end+1)]
    target_nodes_end= target_nodes_end1 + target_nodes_end2
    
    print(target_nodes_start)
    singleSP_start=Cost_radial_start
    singleSP_end=Cost_radial_end
    print(singleSP_start.shape)
    #singleSP_start[1:9,np.arange(89,1,-10)]=1000
    singleSP_start[1:9,:]=1000

    singleSP_start[1:9,10]=0.001

    # singleSP_start[1,89]=1000

    #print(singleSP_start)

    highvalue=int(10000)
    sPath_terminal_start=[]
    sPath_terminal_end=[]   

    #print(target_node_start1,B_path1)

    for i in range(len(target_nodes_start)):

        results=bellman_ford(singleSP_start,N_start,start_node_start,int(target_nodes_start[i]),highvalue)
        sPath_terminal_start.append(results)
       

    for i in range(len(target_nodes_end)):
        sPath_terminal_end.append(bellman_ford(singleSP_end,N_end,start_node_end,target_nodes_end[i],highvalue))

    radial_nodes_long_end=radial_nodes_long_end.tolist()
    radial_nodes_long_start=radial_nodes_long_start.tolist()
    radial_nodes_lat_end=radial_nodes_lat_end.tolist()
    radial_nodes_lat_start=radial_nodes_lat_start.tolist()
    #sPath_terminal_start=sPath_terminal_start.tolist()
    t2=time.time()
    print(t2-t1)

    return JsonResponse ({'radial_nodes_long_start' : radial_nodes_long_start,\
            'radial_nodes_lat_start': radial_nodes_lat_start,'radial_nodes_long_end' : radial_nodes_long_end,\
            'radial_nodes_lat_end': radial_nodes_lat_end,'n_circum_start':n_circum_start,\
            'n_radial_start':n_radial_start,'n_circum_end':n_circum_end,'n_radial_end':n_radial_end,'sPath_terminal_start':sPath_terminal_start,'sPath_terminal_end':sPath_terminal_end})


def shortest_path1(request):

# request is a REST api framework using which we can request variables saved as session in other views modules...
    t1=time.time()
    Cost_horizontal=request.session['Cost_horizontal']
    Cost_vertical=request.session['Cost_vertical']
    Cost_diag=request.session['Cost_diag']
    start=request.session['start']
    Patchx=request.session['Patchx']
    Patchy=request.session['Patchy']
    K=request.session['K']
    L=request.session['L']
    d=request.session['d']
    Lat=request.session['Lat']
    Long=request.session['Long']
    n=int(1/48*((K+2)*(K+4)*(2*K+6)))
    #Node co-ordinates or hashes (we can use any hash generating function here) 
    Hash=np.zeros((n,3))
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
    
    # Dijkstra's cython function call
    Cost_horizontal=np.array(Cost_horizontal)
    Cost_vertical=np.array(Cost_vertical)
    Cost_diag=np.array(Cost_diag)
    Node_val=node_val(Hash,Cost_horizontal,Cost_vertical,Cost_diag,N_K,K,d)

    # Shortest path computation           
    sPath=np.zeros((K,3))    
    sPath[K-1,:]=Node_val[n-1,:]
    for m in range(K-2,-1,-1):
        sPath[m,:]=Node_val[int(sPath[m+1,1]),:]
    #print(sPath)
    
    sPath=sPath.tolist()
    Hash=Hash.tolist()
    #Patchx=Patchx.tolist()
    #Patchy=Patchy.tolist()
    #Lat=Lat.tolist()
    #Long=Long.tolist()
    request.session['N_K']=N_K.tolist()
    request.session['Hash'] = Hash
    
    t2=time.time()
    print(t2-t1)


    return JsonResponse({'sPath':sPath,'start': start,'Hash':Hash,'Patchx':Patchx,'Patchy':Patchy,'Lat':Lat,'Long':Long,'K':K,'d':d})



def Modify_cost(request):
    t1=time.time()

    Patchx=request.session['Patchx']
    Patchy=request.session['Patchy']
    N_K=request.session['N_K']
    N_K=np.array(N_K)  
    N=request.session['N']    
    K=request.session['K']
    d=request.session['d']
    start=request.session['start']  
    Hash=request.session['Hash']      
    Hash=np.array(Hash) 
    
    n=int(1/48*((K+2)*(K+4)*(2*K+6)))
    
    
    Cost_horizontal=request.session['Cost_horizontal']
    Cost_horizontal=np.array(Cost_horizontal)    
    Cost_vertical= request.session['Cost_vertical']
    Cost_vertical=np.array(Cost_vertical)
    Cost_diag= request.session['Cost_diag']
    Cost_diag=np.array(Cost_diag)    
    
    iD=int(request.POST['id_nodes'])
    
    Cost_horizontal[int(iD//N),int(iD%N)+1]=float(request.POST['Cost_horizontal'])
    
    Cost_vertical[int(iD//N)+1,int(iD%N)]=float(request.POST['Cost_vertical'])
    Cost_diag[int(iD//N)+1,int(iD%N)+1]=float(request.POST['Cost_diag'])
    
    sPath1=np.zeros((K,3))
    Node_val=node_val(Hash,Cost_horizontal,Cost_vertical,Cost_diag,N_K,K,d)
    sPath1[K-1,:]=Node_val[n-1,:]
    
    for m in range(K-2,-1,-1):
        sPath1[m,:]=Node_val[int(sPath1[m+1,1]),:]
 
    sPath1=sPath1.tolist()
    Hash=Hash.tolist()
    
    N_K=N_K.tolist()
    Cost_horizontal=Cost_horizontal.tolist()
    Cost_vertical=Cost_vertical.tolist()
    Cost_diag=Cost_diag.tolist()
    

    request.session['Patchx']=Patchx
    request.session['Patchy']=Patchy    
    request.session['N_K']=N_K
    
    request.session['N']=N  
    request.session['K']=K
    request.session['d']=d
    request.session['start']=start
    request.session['Hash']=Hash
      
    request.session['Cost_horizontal'] = Cost_horizontal
    request.session['Cost_vertical'] = Cost_vertical
    request.session['Cost_diag'] = Cost_diag

    t2=time.time()
    print(t2-t1)
     
    return JsonResponse({'sPath1':sPath1,'start': start,'K':K,'d':d,'Hash':Hash,'Patchx':Patchx,'Patchy':Patchy})



