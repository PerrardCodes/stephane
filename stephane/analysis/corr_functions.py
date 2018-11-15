

import numpy as np
import math


import stephane.manager.access as access

import stephane.analysis.statP as statP
import stephane.analysis.cdata as cdata

import stephane.manager.access as access



###################### Correlation functions in space #####################

def corr_d(M,frame,indices=None,dlist=None,axes=['Ux','Ux'],p=1,average=False):
    """
    Compute the correlation function in time at a given instant
    INPUT
    -----
    M : Mdata object
    frame : int
        frame index
    indices : dict of 2 elements tuple (for both keys and values) 
        pairs of coordinates that defines the distance of computation. 
        Default value is None : compute the pair of indices directly
    dlist : list
        distances between points. defaut is None (goes with indices)
    axes : 2 element string list
        field names to be used
    p : int
        order of the correlation function
    OUTPUT
    -----
    dlist : list of distances beetween points
    C : correlation function (un-normalized)
    """
    if indices is None:
        dlist,indices = get_indices(M)
#    else:
#        print('dlist ?')

    C=[]
    X,Y=access.chose_axe(M,frame,axes)
        
    if average:
        Xmoy,Xstd=statP.average(X)
        Ymoy,Ystd=statP.average(Y)
    else:
        Xmoy = 0
        Xstd = 0
        Ymoy = 0
        Ystd = 0

    C = [[] for i in indices]
    for m,ind in enumerate(indices):
        for i,j in ind.keys():  
            k,l= ind[i,j]
            Sp = (X[i,j]-Xmoy)**p*(Y[k,l]-Ymoy)**p   #remove the average in space ?   -> remove by default
            C[m].append(Sp)# for k,l in indices[(i,j)]])    

        C[m],std = statP.average(C[m])
    return dlist,C
    
def corr_d_stat(M,N,indices=None,axes=['Ux','Ux'],p=1,average=False):
    """
    Compute the correlation function from an ensemble N of time
    """
    #select N times to compute the 
    if indices is None:
        dlist,indices = get_indices(M)
        
    nx,ny,nt = M.shape()
    frames = range(0,nt,nt/N)
    
    Ctot = [[] for i in frames]
    for i,frame in enumerate(frames):
         d,Ctot[i] = corr_d(M,frame,indices=indices,dlist=dlist,axes=axes,p=p,average=average)
    
    Ctot = np.asarray(Ctot)
    Ctot,Cstd = statP.average(Ctot,axis=0)
    return dlist,Ctot
    
def get_indices(M,N=10**3):
    """
    Compute N pairs of indices for increasing distances between points
    """
    nx,ny,nt = M.shape()
    dlist=range(int(max([nx/2.,ny/2.])))
    indices=[[] for i in dlist]

    for i,d in enumerate(dlist):
        indices[i]=d_2pts_rand(M.Ux[:,:,0],d,N)   
  #  print('Pair of indices computed')

    return dlist,indices
    
def d_2pts_rand(u,d,N,epsilon=0.5):
    """
    Return N pairs of indice points distant from d  
    """
    #find N pairs of indices such that the distance between two points is d +/- epsilon
    N_rand=10*N
    #from a centered index, return the list of all the matching index
    ny,nx=u.shape
    indices=[]

    i0 = ny//2
    j0 = nx//2
    
    tata,indices_ref=cdata.neighboors(u,(i0,j0),b=d+epsilon/2,bmin=d-epsilon/2)

    i1 = np.floor(np.random.rand(N_rand)*ny).astype(int)
    j1 = np.floor(np.random.rand(N_rand)*nx).astype(int)
    
    theta = np.random.rand(N_rand)*2*math.pi
    i2 = np.floor(i1 + d*np.cos(theta)).astype(int)
    j2 = np.floor(j1 + d*np.sin(theta)).astype(int)
    
    i_in = np.logical_and(i2>=0,i2<ny)
    j_in = np.logical_and(j2>=0,j2<nx)

    keep = np.logical_and(i_in,j_in,np.logical_and(i1==i2,j1==j2))
    
    i1 = i1[keep]
    j1 = j1[keep]
    i2 = i2[keep]
    j2 = j2[keep]
    
    i1 = i1[:N]# keep only the N first indexes
    j1 = j1[:N]
    i2 = i2[:N]
    j2 = j2[:N]

#    print("Number of kept indice pairs : "+str(i1.shape[0]))
    indices = {(i,j):(k,l) for i,j,k,l in zip(i1,j1,i2,j2)}
    
    return indices
    
    
####################### Correlation functions in time ##############################
    
    
