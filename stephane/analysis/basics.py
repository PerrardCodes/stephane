# -*- coding: utf-8 -*-
"""
Created on Wed Oct 21 15:37:51 2015

@author: stephane
"""



import numpy as np
#from scipy import ndimage,ndarray,interpolate
from matplotlib import backend_bases as bck
from math import *
from scipy import ndimage
import os.path
import random
import pylab as plt
import time

import stephane.display.graphes as graphes

def profile_average(M,i,Dt=10,p=1,direction='h'):
    """
    Compute a one  dimensional spatial average of 2d(more ?) data on a time window, centered on a given time.
    INPUT
    -----
    M : Mdata object
        M contains at least Ux and Uy as attributes, and shape() method
    i : int, index of time
        The average will be computed around the i-th time step in M
    Dt : int, time windows on which the average will be computed. Default value : 10
    p : int (or float): order for the p-norm averaging
        p = 1 corresponds to mean
        p = 2 corresponds to standard deviation
    direction : string. Determines the direction of averaging. default value : 'h'
        Can be either 'v' or 'h'. Default value corresponds to a vertical profile as output
    OUTPUT
    -----
    z : 1-dimensional numpy array
        axis values of the 1d profile
    Upx : 1-dimensional numpy array
        averaged field along the direction x      
    Upy : 1-dimensional numpy array
        averaged field along the direction y      
    """
    #start and end corresponds to the bound for the average in time
    #average along the dimensions contained in dimensions for instance : ['x','y','z','t'],
    #the return function is given by the other dimensions (?)
    nx,ny,nt=M.shape()

    start = max(i-Dt//2,0)
    end = min(i+Dt//2,nt-1)
    n=end-start    
    
    Ux=M.Ux[:,:,start:end]
    Uy=M.Uy[:,:,start:end]
    
    if M.param.angle==90:
        axes=(0,2,1) #permute the time and the dimensions whose shouldn't be averaged along (here the x direction in second )
        shape=(nx,n,ny)
        z=M.x[0,:]
    if M.param.angle==0:
      #  print('horizontal average')
        axes=(1,2,0) #permute the time and the dimensions whose shouldn't be averaged along (here the x directionm in second )
        shape=(ny,n,nx)
        z=M.y[:,0]

    Ux=np.reshape(np.transpose(Ux,axes),(shape[0]*shape[1],shape[2]))
    Uy=np.reshape(np.transpose(Uy,axes),(shape[0]*shape[1],shape[2]))

    #RMS velocity profile along the second direction : should correspond to the z axis            
    Ux_moy=np.nanmean(Ux,axis=0)
    Uy_moy=np.nanmean(Uy,axis=0)
    if p==1:
        return z,Ux_moy,Uy_moy
    else:
        Upx=(np.nanmean((Ux-Ux_moy)**p,axis=0))**(1./p)
        Upy=(np.nanmean((Uy-Uy_moy)**p,axis=0))**(1./p)
        return z,Upx,Upy
    
def smooth(V,Dt=10,fixed_length=False):
    """
    Smooth a N+1 dimensional array using a linear average
    The smoothing operation is performed along the lqst dimension (may correspond to time axis)
    INPUT 
    -----
    V : N+1 numpy array to smooth
    Dt : int, number of points for the averaging. Default value is 10
    fixed_length : not implemented yet
    
    OUTPUT
    -----
    ts : 1 dimensional array,
        new axis for the dimension of the smoothing
    Vs : smoothed N+1 dimensional array
        The length along the averaged dimension has been shortened by Dt
    """
    dim=np.shape(V)
    print(dim)
    nt=dim[-1]    
    N=len(dim)   
    
    if fixed_length:
        print("add zero matrices at the begining and the end to keep the length constant")
        print("Not implemented yet")
        
  #  print(Dt)
  #  print([k for k in range(Dt,nt-Dt)])
    Vs=np.array([np.nanmean(V[...,k-Dt:k+Dt],N-1) for k in range(Dt,nt-Dt)])
        
    #some elements are lost during the smoothing operation : we should add the initial and enn elements
 #   Vs=np.array([V[...,:Dt]+[np.nanmean(V[...,k-Dt:k+Dt],N-1) for k in range(Dt,nt-Dt)])

#    tup=
    order=np.arange(1,N+1)
    order[N-1]=0
    ind=tuple(order)
  #  print(Vs.shape)
  #  print(ind)
    
    Vs=np.transpose(Vs,ind)
#    indices=[slice(dim[i]) for i in range(N-1)]  
  #   indices[N-1]=
#    for k in range(Dt,nt-Dt):
 #       Vs.append())

#    if len(dim)==1:
#        Vs=np.array([np.mean(V[k-Dt:k+Dt]) for k in range(Dt,nt-Dt)])
#    if len(dim)==2:
#       Vs=np.array([np.mean(V[:,k-Dt:k+Dt],N-1) for k in range(Dt,nt-Dt)])
#    if len(dim)==3:
#        Vs=np.array([np.mean(V[:,:,k-Dt:k+Dt],2) for k in range(Dt,nt-Dt)])
#    if len(dim)==4:
#        Vs=np.array([np.mean(V[:,:,:,k-Dt:k+Dt],3) for k in range(Dt,nt-Dt)])
    return Vs

def fluctuations(S,Dx,Dy,Dt):
    """  
    compute the standard deviation from a volume (Dx,Dy,Dt) around each point.
    Return an array of same size than Vec, containing the local standard deviation
    
    """
    U=S.Ux**2+S.Uy**2
    Ufilt2=ndimage.filters.gaussian_filter(U,[Dx,Dy,Dt])
   
    Z=U-Ufilt2
    start=0
    stop=1500
    Dirname='Velocity_fluctations'
    make_2dmovie(S,Z,start,stop,Dirname)
    
def horizontal_profile(S,ylines,Dt,start=0):
    nx,ny,nt=S.shape()
    
    x=S.x[0,:]
    for i in range(start,nt,Dt):
        Ux=np.mean(np.mean(S.Ux[ylines,:,i:i+Dt],axis=0),axis=1)
        Uy=np.mean(np.mean(S.Uy[ylines,:,i:i+Dt],axis=0),axis=1)
        
        std_Ux=np.std(np.std(S.Ux[ylines,:,i:i+Dt],axis=0),axis=1)
        std_Uy=np.std(np.std(S.Uy[ylines,:,i:i+Dt],axis=0),axis=1)
        
        plt.subplot(121)
        graphes.graph(x,Ux,0,std_Ux)
        graphes.legende('x (m)','V (m/s)','Ux')

        plt.subplot(122)
        graphes.graph(x,Uy,0,std_Uy)
        graphes.legende('x (m)','V (m/s)','Uy')
        
        plt.draw()
        raw_input()

def vertical_profile(S,xlines,Dt,start=0):
    nx,ny,nt=S.shape()
    
    y=S.y[:,0]
    for i in range(start,nt,Dt):
        Ux=np.mean(np.mean(S.Ux[:,xlines,i:i+Dt],axis=1),axis=1)
        Uy=np.mean(np.mean(S.Uy[:,xlines,i:i+Dt],axis=1),axis=1)
        
        #standard deviation computation
        std_Ux=np.sqrt(np.mean(np.mean(abs(S.Ux[:,xlines,i:i+Dt]-Ux)**2,axis=1),axis=1))
        std_Uy=np.sqrt(np.mean(np.mean(abs(S.Uy[:,xlines,i:i+Dt]-Uy)**2,axis=1),axis=1))
        
        print(std_Ux)
        
        plt.subplot(121)
        graphes.graph(y,Ux,std_Ux)
        graphes.legende('z (m)','V (m/s)','Ux')

        plt.subplot(122)
        graphes.graph(y,Uy,std_Uy)
        graphes.legende('z (m)','V (m/s)','Uy')
        
        plt.draw()
        raw_input()
    
def compare_profil(S1,S2):
    #### DEPRECIATED : move Sdata.subset_index to a new module of smart data indexing

    #S1 and S2 must be the same length, and the same x and y dimensions
    #compare velocity measurement obtained with two different frame rate
    #norm and direction
    Ux1,Uy1=fix_PIV(S1)
    Ux2,Uy2=fix_PIV(S2)
    
    U1,theta1=Smath.cart2pol(Ux1,Uy1)
    U2,theta2=Smath.cart2pol(Ux2,Uy2)
    print(U2.shape)
    
    nx,ny,nt=S1.shape()
    
    indices1,indices2,nt=subset_index(S1,S2)
    #locate a subset of identical im_index
    #indices correspond of the number of the U rows, not to the im_index valeurs
  #  indices,nt=subset_index(S1,S2)
    
    U1=np.reshape(U1[:,:,indices1],nx*ny*nt)
    U2=np.reshape(U2[:,:,indices2],nx*ny*nt)
    
    #if U1 and U2 have not the same length,
    #a fraction can be selected using the attr im_index of each one (to use the same list of images)
    #random extract
    N=2000
    ind=random.sample(xrange(nx*ny*nt),N)
    
    graphes.graph(U1[ind],U2[ind])
    xlabel='V (m/s)  at '+str(S1.timescale*S1.fps)+' fps'
    ylabel='V (m/s)  at '+str(S2.timescale*S2.fps)+' fps'

    graphes.legende(xlabel,ylabel,'')

    bounded_velocity(S1,True,5,'v')
    bounded_velocity(S2,True,5,'h')
    
    raw_input()    
    
def subset_index(S1,S2):
    return None

def mean_profile(S,i,j,direction='v',label='k^',display=False):
    #mean profile along the whole field : average on one direction only ! (and small windows on the other direction ?)
    nx,ny,nt=S.shape()

    Ux=S.m.Ux
    Uy=S.m.Uy
    #remove the data out of the PIV bounds    
#    Ux,Uy=fix_PIV(S)

    U=np.sqrt(Ux**2+Uy**2)
#    V=np.reshape(U,(nx*ny,nt))
    #median is not so affected by peak values, but standard deviation definetely !
    #histogramm between vmin and vmax, and remove values out of bound (set to NaN)
    U_moy=[]
    U_std=[]
    t=S.m.t

    Dt=2
    if direction=='v':
        #average along the horizontal direction
        U_moy=[np.mean(np.mean(U[j-Dt:j+Dt,:,k],axis=0),axis=0) for k in range(nt)]
        print('horizontal average')
    else:
        #average along the vertical direction
        U_moy=[np.mean(np.mean(U[:,i-Dt:i+Dt,k],axis=0),axis=0) for k in range(nt)]
        print('vertical average')
        
    print(np.shape(U_moy))
    if display:
    #U_moy=np.mean(V[np.invert(np.isnan(V))],axis=0)        
        print('Number of frames : '+str(len(S.m.t)))
    
        graphes.graph(t,U_moy,label)
        graphes.legende('t (ms)','<V>_{x,y} (m/s)','')
    
    return U_moy,U_std    
#    plt.xlim([0, 500])
#    plt.ylim([0, 3])
    

def profile_xy(S,x0,y0):
    #plot the profile at a given position : will be use to compare several movies
    i,j=get_index(S.m,x0,y0)
    return i,j
    
def histogramm(M,label='ko'):
    #complete histogramm
    nx,ny,nt=M.shape()
    U=np.sqrt(M.Ux**2+M.Uy**2)
    
    V=np.arange(0,2,0.005)
    n,x,shape=plt.hist(np.reshape(U,nx*ny*nt),V,histtype='step')
    
    #add bounded values of velocity measurement ??
#    bounded_velocity(M,True,max(n))
  #  plt.show(block=False)
  #  time.sleep(10e-3)
        

def velocity_fluctuation(S,nframe):
    #for each point, compute an everage variation of velocity on nframe
    #both in angle and modulus !
    nx,ny,nt=S.shape()

    dUx=np.diff(S.Ux)
    dUy=np.diff(S.Uy)
    
    dU=dUx**2+dUy**2
    
    #cumulative sum on the increment of velocity
    dU_sum=np.cumsum(dU,axis=2)
 
    dU_nframe=np.sqrt((dU_sum[:,:,nframe:]-dU_sum[:,:,0:-nframe])/nframe)
    
    make_2dmovie(S,dU_nframe,1,nt-nframe,'Velocity_increment')
    
#    dU_moy=[dU_sum[:,:,i]-dU_sum[:,:,0] for i in range(1,nframe)]
#    dUr,dUtheta=Smath.cart2pol(dUx,dUy)
    
def velocity(M,nframe=0):
    nx,ny,nt=M.shape()
    
    U=(M.Ux**2+M.Uy**2)
 
    if nframe==0:
        U_nframe=U
    else:    
        U_sum=np.cumsum(U,axis=2)        
        U_nframe=np.sqrt((U_sum[:,:,nframe:]-U_sum[:,:,0:-nframe])/(nframe+1))
    
#    Dir_root = os.path.dirname(M.Sdata.fileCine)
    Dirname='Velocity_modulus_'+M.Id.get_id()
    make_2dmovie(M,U_nframe,1,630,Dirname)


def time_window(S,nframe):
    #compute histogramm from each time windows
    U=np.sqrt(S.Ux**2+S.Uy**2)
    nx,ny,nt=S.shape()
    #average every nframe :
    for i in range(nframe,nt,nframe):
        V=U[i-nframe:i]
        
        plt.clf()
        n,x,shape=plt.hist(np.reshape(V,nx*ny*nframe),V,histtype='step')
        bounded_velocity(S,True,max(n))
            
        plt.show(block=False)
        time.sleep(10e-3)
        raw_input()
    
    
def multi_graph(function,Slist):
    for S in Slist:
        function(S)
        raw_input()

def multi_mean_profile(Slist,i,j,direction='v',shift=0):
    label=['k^','ro','bx','mo','gx','k+']
    for S in Slist:
        plt.clf()
        mean_profile(S,i,j,direction,label[(Slist.index(S)+shift)%len(label)],True)
#        mean_profile(S,label[Slist.index(S)%len(label)])

def fix_PIV(S):
    #replace by NaN all the values of Ux and Uy that does not correspond to a bounded velocity
    Umax,Umin=bounded_velocity(S)
    
    U=np.sqrt(S.Ux**2+S.Uy**2)

    Vtrue=np.logical_and(U>Umin,U<Umax)
 
    Ux=np.where(Vtrue,S.Ux,np.NaN)
    Uy=np.where(Vtrue,S.Uy,np.NaN)
    
    Ntrue=sum(sum(sum(np.where(Vtrue,1,0))))
    Nfalse=sum(sum(sum(np.where(Vtrue,0,1))))
    N=Ntrue+Nfalse
    
    print('Percentage of removed data points : '+str((Nfalse*100)/N)+ ' %')
    
    return Ux,Uy
 
def divergence(S,xlines,ylines,display=False):
    """
    DEPRECIATED : was used to compute the asymetry of the flow
    
    """
    #xlines list of tuple ??
    #take two points (symmetric in respect to X=0?), and mesure the relative distance between each
    #as a function of time
    #relative means in respect to the spatial averaged velocity at the time
    t=S.t
    Dt=10
    ny,nx,nt=S.shape()
    
    U_moy,U_std=mean_profile(S,label='k^',display=False)
    
    n=len(xlines)
    Ds=np.zeros((ny,nx,nt-2*Dt)) #first element is centered !
    
    for i in ylines:
        for pair in xlines:
            dx=S.Ux[i,pair[0],:]+S.Ux[i,pair[1],:] #antisymetric x profile ?
            dy=S.Uy[i,pair[0],:]-S.Uy[i,pair[1],:]
            D,phi=Smath.cart2pol(dx,dy)

            #Divide by the global average velocity at each time !!    
            Ds[ylines.index(i),pair[0],:]=smooth(D,Dt)#=smooth(D/U_moy,Dt)
            Ds[ylines.index(i),pair[1],:]=Ds[i,pair[0],:]
            
            if display:
                graphes.graph(t[Dt:-Dt],Ds)
                graphes.legende('t (ms)','V (m/s)','Velocity difference between two symetric points')
    return Ds


    #not NaN elements can be acceded from 
#    np.median(a[np.invert(np.isnan(a))])
#def main():
 #   velocity(M_log[2])

#main()    
    