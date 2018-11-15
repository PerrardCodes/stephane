import numpy as np
import matplotlib.pyplot as plt
import datetime
import os
import math
import argparse
import sys

import os.path
from matplotlib.patches import Circle, Wedge, Polygon
from matplotlib.collections import PatchCollection
import scipy
from mpl_toolkits.mplot3d import Axes3D

import stephane.hdf5.h5py_s as h5py
import stephane.display.graphes as graphes
    
import time

import stephane.design_vortices.random as rand_tangle

##################### Initialization ##################
        
def pvec(t):
    """
    pick one "arbitrary" perpendicular vector of t
    """
    e1 = [1,0,0]
    e2 = [0,1,0]
    e3 = [0,0,1]
    
    u = np.cross(t,e1)+np.cross(t,e2)#+np.cross(t,e3)
#    print(u)
    if np.linalg.norm(t)==0:
        print("zero vector !")
        
    return u/np.linalg.norm(u)    
    
    
def Rop(t,npoints,d=3):
    """
    Generate a Rotation operator along the plane perpendicular to t
    """
    
    Theta = np.arange(0,2*np.pi,2*np.pi/npoints)
    
    Tcos_mat = np.transpose(np.tile(np.cos(Theta),(3,3,1)),(2,0,1))
    Tsin_mat = np.transpose(np.tile(np.sin(Theta),(3,3,1)),(2,0,1))

    e1 = [1,0,0]
    e2 = [0,1,0]
    e3 = [0,0,1]
    t_cross = np.asarray([np.cross(e1,t),np.cross(e2,t),np.cross(e3,t)])
    
    
    Id = Tcos_mat*np.tile(np.identity(d),(npoints,1,1))
    Cross = Tsin_mat*np.tile(t_cross,(npoints,1,1))
    Tensor = (1-Tcos_mat)*np.tile(np.tensordot(t,t,axes=0),(npoints,1,1))
    
    return Id+Cross+Tensor

################### Tools ###################

def tangent(X,d=3,step=1,cyclic=True):
    """
    Compute the tangent vector on each point. 
        Assume that the first and the last point are close each other (closed loop !)
        how should we normalize it ? -> by the curviligne abscisse 
    INPUT
    -----
    X : [N,d] array
        Line of spatial coordinates
    d : int
        spatial dimension, default is 3
    step : int
        index step used to compute the spatial derivative. default is 1
    OUTPUT
    -----
    dV : 
    """
    alpha=3/4.
    beta=-3/20.
    gamma=1/60.
    
#    alpha = 1.
#    beta = 0.
#    gamma = 0.

    N = np.shape(X)[0]
    n=3    
    
    if cyclic:
        X_ext = np.concatenate((X[N-n:,...],X,X[0:n,...]),axis=0)
        dV = np.zeros((N,d))
    else:
        X_ext = X
        dV = np.zeros((N-n*2,d))
        
    if d==1:
        X_ext = np.reshape(X_ext,(np.shape(X_ext)[0],1))
    
    #compute the first derivative along each spatial axis
    for j in range(d):
        
        tl=[[slice(3,-3,step)]+[j] for p in range(6)]
    
        tl[0][0]=slice(0,-5,step)
        tl[1][0]=slice(1,-4,step)
        tl[2][0]=slice(2,-3,step)
        tl[3][0]=slice(3,-2,step)
        tl[4][0]=slice(4,-1,step)
        tl[5][0]=slice(5,None,step)
    
        dV1 = np.sum(np.asarray([np.diff(X_ext[tl[k]],axis=0) for k in range(2,4)]),axis=0)
        dV2 = np.sum(np.asarray([np.diff(X_ext[tl[k]],axis=0) for k in range(1,5)]),axis=0)
        dV3 = np.sum(np.asarray([np.diff(X_ext[tl[k]],axis=0) for k in range(6)]),axis=0)
        
        dV[...,j] = alpha*dV1+beta*dV2+gamma*dV3
    return dV


def norm(u,axis=1):
    return np.sqrt(np.sum(u**2,axis=axis))


def normalize(U):
    """
    Normalize the matrix U along its last dimension
    """
    d = len(np.shape(U))-1
    
    norm = np.sqrt(np.sum(U**2,axis=d))
    permute = tuple([i for i in range(1,d+1)]+[0])
    
    U_norm = np.transpose(np.asarray([norm for k in range(d+1)]),permute)
    
    return U/U_norm
    
########################### Start ###################

def generate_rings(Radius,positions,vectors,npoints):
    """
    generate a serie of rings at the specified locations positions with an oriented vector  vectors
    """
    
    theta = np.arange(0,2*np.pi,2*np.pi/npoints)
        
    U = []
    for pos,vector in zip(positions,vectors):        
        R = Rop(vector,npoints) 
       # print(vector)
        u_n = pvec(vector)
        u = np.tile(pos,(npoints,1)) + np.dot(R,u_n)*Radius   # a could be a function of Theta ...
        U.append(u)
    return U
    
def example(npoints):
    positions = [[0,1,0.5]]#,[0,0,-1],[0,1,0],[0,-1,0]]    
    R = 0.2
    U = generate_rings(R,positions,positions,npoints)
    
    return U
    #    plt.plot(U)
    
def circle(nring,R=1):
    theta = np.arange(0,2*np.pi,2*np.pi/nring)
    positions = [[R*np.cos(t),R*np.sin(t),0] for t in theta]
  #  print(np.asarray(positions).shape)
    
    return positions
    
def add_vector(Dict,R):
    """
    add a vector field to Dict['R']
    """

    for i,Rp in enumerate(Dict['Rp']):
        Dict['Rp'][i]=Dict['Rp'][i]+R[i]

        
def add_noise(Dict,sigma,n=10,recompute=True,An=None):
    t,An = rand_tangle.noise(Dict['U'],sigma,n,T=0,N=Dict['npoints'],remove_mean=False,recompute=recompute,An=An)  
      
  #  print(t.paths[0].shape)
    Dict['U']=t.paths[0]
    Dict['An']=An
    
    return Dict
        
############################ Computation #####################

    
def compute_serie(Dict,Gammalist,d=3):
    # to be fast 
    R = Dict['R']
    m = R.shape[0]

    #memory allocations :
    #C, T_mat, norm_C
    
    epsilon = Dict['epsilon']
    Dict['U'] = np.zeros((m,3))

    for Rp,Gamma in zip(Dict['Rp'],Gammalist):
        n = Rp.shape[0]

        T = tangent(Rp,d=3,step=1,cyclic=True) #compute tangent vector
        Dict['T_mat'] = np.transpose(np.tile(T,(m,1,1)),(1,0,2)) # make a (n,m,3) matrix
        
        Dict['C'] = np.tile(R,(n,1,1))-np.transpose(np.tile(Rp,(m,1,1)),(1,0,2))        
        Cn = norm(Dict['C'],axis=2)
        Dict['norm_C'] = np.transpose(np.asarray([Cn for k in range(d)]),(1,2,0))+epsilon
    
        Dict['U'] += Gamma/(4*np.pi)*np.sum(np.cross(Dict['T_mat'],Dict['C'])/Dict['norm_C']**3,axis=0)
    
    return Dict
    
def initialize(nring,npoints):
    """
    Initializing function. Write all the parameters in a dictionnary.
    Those parameters will be saved with the result of the simulation
    """
    Dict = {}
    
    n = npoints#V.shape[0] # number of points
    m = npoints*nring#U.shape[0] #total number of points

    Dict['dt'] = 0.001
    Dict['radius']=1.
    Dict['Radius']=1.
    Dict['epsilon']=0.1#025
    Dict['noise']=1.#10**-12
    Dict['noise_increment']=0.01
    Dict['n_mode']=10
    Dict['C']= np.zeros((n,m,3))
    Dict['T_mat']= np.zeros((n,m,3))
    Dict['norm_C']= np.zeros((n,m,3))
    
    Dict['nring'] = nring
    Dict['npoints'] = npoints
    
    
  #  Dict['Rp'] = example(npoints)#np.zeros((npoints,3))
    positions = circle(nring,R=Dict['Radius'])#np.zeros((npoints,3))
    Dict['Rp'] = generate_rings(Dict['radius'],positions,positions,npoints)
    
 #3   print(np.asarray(Dict['Rp']).shape)
    
    Dict['R'] = np.reshape(np.asarray(Dict['Rp']),(npoints*nring,3))#np.zeros((npoints,3))
 #   print(Dict['R'].shape)
    
    Dict['U'] = np.zeros((m,3))
#    print(Dict['U'].shape)
#    U['paths'] = generate_rings(R,N,positions)    
    return Dict

        
def compute(N,noise=False,display=False):    
    npoints = 200
    #t1 = time.time()
    nring = 1
    Dict = initialize(nring,npoints)
    sigma = 0.1
    n=5
    t,An = rand_tangle.noise(Dict['U'],sigma,n,T=0,N=Dict['npoints'],remove_mean=False)  
    Dict['R']=Dict['R']+t.paths[0]
    
    Dict = add_noise(Dict,sigma,n=Dict['n_mode'])        
    
    #t2 = time.time()
    #print(t2-t1)
    c=0
    

    keys = ['Radius','radius','dt','epsilon','nring','npoints','noise','n_mode']
    
    s=""
    for key in keys:
        s=s+key+'_'+str(Dict[key])+'_'
    s=s[:-1]    
    folder = s
    
    Dict['R_time']=[None for j in range(nring)]
    for j in range(nring):    
        Dict['R_time'][j] = np.zeros((npoints,3,N))#Dict['R'][j*npoints:(j+1)*npoints]
    
    
    
    for i in range(N):
     #   if i%1000==0:
     #       print(i)
     #   t1 =time.time()
        Dict = compute_serie(Dict,1*np.ones(nring),d=3)
    
        if noise:
            mean = np.sqrt(np.mean(Dict['U']**2))
        #print(mean)
            sigma = mean*Dict['noise']
       # print(sigma)
       # print(sigma
        #if i==0:
            Dict = add_noise(Dict,sigma,n=Dict['n_mode'])        
        #else:
        #    Dict['An']=np.exp(-Dict['noise_increment'])*Dict['An']
        #    Dict = add_noise(Dict,sigma,n=Dict['n_mode'],recompute=False,An=Dict['An'])
        #    Dict = add_noise(Dict,sigma*(1-np.exp(-Dict['noise_increment'])),n=Dict['n_mode'])
         #   if i<10:

        Dict['R'] = Dict['R']+ Dict['dt']*Dict['U']
        
        for j in range(nring):    
            Dict['Rp'][j] = Dict['R'][j*npoints:(j+1)*npoints]
            Dict['R_time'][j][...,i] = Dict['R'][j*npoints:(j+1)*npoints]
        
        if display:
            if i==0:
                fig = plt.figure(1)
                ax = fig.add_subplot(111, projection='3d')
                
                for Rp in Dict['Rp']:
                    ax.plot(Rp[:,0],Rp[:,1],Rp[:,2])
                plt.axis('equal')
                
                #plt.plot(Dict['R_time'][0][:,1],Dict['R_time'][0][:,2])
                #plt.show()
        
        
#        if i%10==1:
            #plt.cla()            

#            name = "N="+str(nring)+" rings"
#            figs = graphes.legende('X','Y',"N="+str(nring)+" rings",cplot=True)
#            graphes.save_figs(figs,savedir='./Numeric/Biot_Savart/'+folder+'/',suffix=str(i),prefix='',frmt='png',dpi=300,display=False)          
    
    return Dict
    #    t2 = time.time()
    #    print(t2-t1)

def length(Dict,i):
    
    R = Dict['R_time'][0][...,i]
    dV = tangent(R,d=3,step=1,cyclic=True)
    
    return np.sum(norm(dV))
    

def main():
    N = 10**7
    num_gyro = 16
    simulate(N,num_gyro)
    
    
if __name__ == '__main__':
    main()
    
    
    ################## Garbage ####################
    

def compute_u_matrix(R,Rp,Gamma,d=3):
    
    #to understand what is going on
    epsilon = 0.05
    
    n = Rp.shape[0]
    m = R.shape[0]

    C = np.tile(R,(n,1,1))-np.transpose(np.tile(Rp,(m,1,1)),(1,0,2))

    T = biot.tangent(Rp,d=3,step=1,cyclic=True) #compute tangent vector
    T_mat = np.transpose(np.tile(T,(m,1,1)),(1,0,2)) # make a (n,m,3) matrix
    
    dl = np.cross(T_mat,C)
    
    Cn = norm(C,axis=2)
    norm_C = np.transpose(np.asarray([Cn for k in range(d)]),(1,2,0))+epsilon
    
    return Gamma/(4*np.pi)*np.sum(dl/norm_C**3,axis=0)
    
def compute_fast(Dict,Gamma,d=3):
    # to be fast 
    R = Dict['R']
    Rp = Dict['Rp']
    #memory allocations :
    #C, T_mat, norm_C
    
    epsilon = 0.1
    n = Rp.shape[0]
    m = R.shape[0]

    Dict['C'] = np.tile(R,(n,1,1))-np.transpose(np.tile(Rp,(m,1,1)),(1,0,2))

    T = tangent(Rp,d=3,step=1,cyclic=True) #compute tangent vector
    
    Dict['T_mat'] = np.transpose(np.tile(T,(m,1,1)),(1,0,2)) # make a (n,m,3) matrix
        
    Cn = norm(Dict['C'],axis=2)
    Dict['norm_C'] = np.transpose(np.asarray([Cn for k in range(d)]),(1,2,0))+epsilon
    
    Dict['U'] = Gamma/(4*np.pi)*np.sum(np.cross(Dict['T_mat'],Dict['C'])/Dict['norm_C']**3,axis=0)
    
    return Dict    
