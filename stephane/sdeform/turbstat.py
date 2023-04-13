import matplotlib.pyplot as plt
import glob
import numpy as np
import stephane.display.graphes as graphes
import pickle
import os

from scipy.spatial import SphericalVoronoi

from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from mpl_toolkits.mplot3d import proj3d

from matplotlib import colors
import scipy

import stephane.sdeform.spherical as sphere


def load_data(files,N=128):
    
    d={}

    for i,filename in enumerate(files):
        dat = np.fromfile(filename, dtype=float)

        nel = len(dat)
        m = int(nel/(N+1)**2)
#print(m)
        dat = np.reshape(dat,(m, N+1, N+1))
    
        d['x'] = dat[2:,0,0]
        d['z'] = dat[0, 1:, 0]
        d['y'] = dat[0, 0, 1:]
        key = os.path.basename(filename)[:2]
        d[key] = dat[2:, 1:, 1:]
    #print(d.keys())

    velocity = ['ux','uy','uz']

    for key in velocity:
        d[key] = d[key] - np.mean(d[key],axis=(0,1,2))
    
    fig = graphes.hist(np.reshape(d['ux'],N**3),log=True)
    fig = graphes.hist(np.reshape(d['uy'],N**3),log=True)
    fig = graphes.hist(np.reshape(d['uz'],N**3),log=True)

    ustd = np.std(np.reshape(d['ux'],N**3))

    uth = np.linspace(-20,20,10**3)
    plt.semilogy(uth,1/np.sqrt(2*np.pi)/ustd*np.exp(-(uth/ustd)**2/2),'r--')
    
    return d

def load_pressure(filename,d,N=128):
    dat = np.fromfile(filename, dtype=float)

    nel = len(dat)
    #print(nel)
    m = int(nel/(N+1)**2)
#print(m)
    dat = np.reshape(dat,(m, N+1, N+1))
    d['p'] = dat[2:, 1:, 1:]
    #print(d.keys())
    
    return d
    
def corr(D,p=2,d=3,N=128):
    # return only the second moment
    Nd = int(N/2)
    
    dU = np.zeros((Nd,3))
    
    for i in range(d):
        for b in range(1,Nd):#b=0 is by definition 0
        #print(b)            
            tuple1 = (slice(b,N,1),slice(None,None,1),slice(None,None,1))
            tuple2 = (slice(0,N-b,1),slice(None,None,1),slice(None,None,1))

            tup1 = tuple(np.roll(tuple1,i))
            tup2 = tuple(np.roll(tuple2,i))
            
            dU[b,i] = np.mean((D[tup1]-D[tup2])**p,axis=(0,1,2))
            
    return dU
    
def corr_stat(D,b=8,d=3,p=2,N=128):
    #return the full distribution
    dU = []
    
    for i in range(d):
        tuple1 = (slice(b,N,1),slice(None,None,1),slice(None,None,1))
        tuple2 = (slice(0,N-b,1),slice(None,None,1),slice(None,None,1))

        tup1 = tuple(np.roll(tuple1,i))
        tup2 = tuple(np.roll(tuple2,i))
        
        (Nx,Ny,Nz) = D[tup1].shape
        dU.append(np.squeeze(np.reshape((D[tup1]-D[tup2])**p,(Nx*Ny*Nz,1,1))))

    #N = len(dU[0])
    #dU = np.reshape(dU,N*3)
    return dU
    
def spherical_increments(M,t,R0,R,Theta=[],Phi=[],A=[],No=30):
    (x0,y0,z0) = R0
    if Theta==[]:
        Theta,Phi = sphere.mesh(No)
    Nv = np.asarray([np.sin(Theta)*np.cos(Phi),np.sin(Theta)*np.sin(Phi),np.cos(Theta)])

    if A==[]:
        [p,sv,A] = sphere.shape(0,0,0,R0=R0,No=No,harm0=R)
        
    X = R*np.sin(Theta)*np.cos(Phi)+x0
    Y = R*np.sin(Theta)*np.sin(Phi)+y0
    Z = R*np.cos(Theta)+z0

    Np = Theta.shape[0]

    U = M['f'](t,(X,Y,Z))
    Umoy = np.sum(A*U,axis=1)/4/np.pi
    Utile = np.transpose(np.tile(Umoy,(Np,1)))
    dU_sphere = np.sum(Nv*(M['f'](t,(X,Y,Z))-Utile),axis=0) 
    
    if 'fp' in M.keys():
        P = M['fp'](t,(X,Y,Z))
        Pmoy = np.sum(A*P)/4/np.pi    
        dP_sphere = np.sum(Nv*(M['fp'](t,(X,Y,Z))-Pmoy),axis=0) 
    else:
        dP_sphere=None
        
    return dU_sphere,dP_sphere

def moment_sphere(U,A,p=2):
    last = len(U.shape)-1
    return np.mean((np.sum(A*U**2,axis=last)/4/np.pi),axis=range(last-1))**(1/2)
