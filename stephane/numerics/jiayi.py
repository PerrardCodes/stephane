import numpy as np
import matplotlib.pyplot as plt
import datetime
import os
import plotmovie as pm
import math
import matplotlib
import argparse
import sys
import cPickle as pickle
import motor_track_functions as mtf
import os.path
from matplotlib.patches import Circle, Wedge, Polygon
from matplotlib.collections import PatchCollection
import scipy
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import jiayi.hdf5.h5py_convert as h5py
    

def laplacian(u,j,dx):
    n=len(u[:,j])
    ddu = (-2*u[1:-1,j]+u[:-2,j]+u[2:,j])/dx**2
    
    ddu_list = np.ndarray.tolist(ddu)
    ddu_list = [ddu[-1]]+ddu_list+[ddu[0]]
    
    ddu = np.asarray(ddu_list)
    
    return ddu

def laplacian_E(u,dx):
    ddu = (-2*u[1:-1]+u[:-2]+u[2:])/dx**2
    return ddu

def iterate_time(U,j):
    """
    U is a dictionnary with fields :
        uij
        phij
        vuij
        vphij
        
        p : dictionnary with fields :
            wspin
            c
            w
            spatialstep
            timestep
            time
            num_gyros
    """
    secondDy = laplacian(U['uij'],j,U['p']['spatialstep'])
    secondDx = laplacian(U['phij'],j,U['p']['spatialstep'])
    
    inertiaphi = U['p']['w']**2*U['phij'][:,j]
    inertiar = U['p']['w']**2*U['uij'][:,j]
    
    #### Forcing term in time
    omega= U['p']['w']
    A=0.
    
    n = U['uij'].shape[0]
    Fu = 0.#A*np.sin(omega*t)
    t = j*U['p']['timestep']
    #Fphi = 0.
    Fphi = np.zeros(n)
    Fphi[4]= np.mean(A*np.cos(np.linspace(U['p']['w'],U['p']['w']*10,100)*t))
    
    ##nonlinearity
    a=5.
    #print(Fphi)
    #print(U['vuij'][:,j+1] + Fphi)
    
    #### Update position from previous velocity 
    U['uij'][:,j+1] = U['uij'][:,j] + U['p']['timestep']*U['vuij'][:,j]
    U['phij'][:,j+1] = U['phij'][:,j] + U['p']['timestep']*U['vphij'][:,j]
    
    #### Update velocity from the equation of motion
    U['vuij'][:,j+1] = U['vuij'][:,j] + U['p']['timestep']*(-secondDy*U['p']['c']/4 - inertiar + U['p']['wspin']*U['vphij'][:,j]+Fu)
    U['vphij'][:,j+1] = U['vphij'][:,j] + U['p']['timestep']*(secondDx*U['p']['c'] - inertiaphi - U['p']['wspin']*U['vuij'][:,j]+Fphi+a*U['phij'][:,j]**3)
    
    #print(U['vphij'][:,j+1])
    ### take care of boundaries
    """
    U['uij'][0,j+1]=U['uij'][-2,j+1]
    U['uij'][-1,j+1]=U['uij'][1,j+1]
    U['phij'][0,j+1]=U['phij'][-2,j+1]
    U['phij'][-1,j+1]=U['phij'][1,j+1]
    U['vuij'][0,j+1]=U['vuij'][-2,j+1]
    U['vuij'][-1,j+1]=U['vuij'][1,j+1]
    U['vphij'][0,j+1]=U['vphij'][-2,j+1]
    U['vphij'][-1,j+1]=U['vphij'][1,j+1]
    """
    #print(U['uij'][0,j+1]-U['uij'][-2,j+1])
    
    ### Normalize total energy
    E = energy(U,j)
    E_new = energy(U,j+1)
    
    beta = np.sqrt(E/E_new)
    
    U['uij'][:,j+1] = U['uij'][:,j+1]*beta
    U['phij'][:,j+1] = U['phij'][:,j+1]*beta
    U['vuij'][:,j+1] = U['vuij'][:,j+1]*beta
    U['vphij'][:,j+1] = U['vphij'][:,j+1]*beta
    
#    A_u = np.sqrt(np.sum(U['uij'][:,j]**2))
#    A_phi = np.sqrt(np.sum(U['phij'][:,j]**2))
    
#    A_u_new = np.sqrt(np.sum(U['uij'][:,j+1]**2))
#    A_phi_new = np.sqrt(np.sum(U['phij'][:,j+1]**2))
    #print(U['uij'])
    return U


def energy(U,j):
    #Energy only from the spring : should add the laplacian in space
    Ep = U['p']['w']**2*np.sum(U['uij'][1:-1,j]**2+U['phij'][1:-1,j]**2)
    
    E_c1 = 1/4*np.sum(laplacian_E(U['uij'][:,j]**2,U['p']['spatialstep']))
    E_c2 = -np.sum(laplacian_E(U['phij'][:,j]**2,U['p']['spatialstep']))
    E_coupling = U['p']['c']*(E_c1+E_c2)
    
    #kinetic energy with m=1 !
    Ec = np.sum(U['vuij'][1:-1,j]**2+U['vphij'][1:-1,j]**2)
    Etot = Ec + Ep + E_coupling
    
    return Etot
    
def initialize(num_gyros,time):
    
    U = {}
    U['uij'] = np.zeros((num_gyros+2,time),dtype=np.float128)
    U['phij'] = np.zeros((num_gyros+2,time),dtype=np.float128)
    U['vuij'] = np.zeros((num_gyros+2,time),dtype=np.float128)
    U['vphij'] = np.zeros((num_gyros+2,time),dtype=np.float128)
    
    U['phij'][num_gyros/2,0]= 0.5
#    U['vuij'][2,0]=0.1
    
    g = 9.81
    l = 0.038
    
    U['p']={}
    U['p']['wspin'] = 0.
    U['p']['c'] = 15 #interaction.
    U['p']['w'] = 10.#np.sqrt(g/l)
    U['p']['spatialstep'] = 1.#0.75
    U['p']['timestep'] = 0.001#25
    U['p']['num_gyros'] = num_gyros
    U['p']['time'] = time
    
    return U

def simulate(N,num_gyros):
#    time = 10
    
    U = initialize(num_gyros,N)
    
    for j in range(N-1):
        if j%(N/10)==0:
            print(j)
        U = iterate_time(U,j)

    num_gyros=U['p']['num_gyros']
    current_time = datetime.datetime.now().strftime("%B_%d_%Y_%I_%M_%p")
    dat_hd = '/Volumes/Jiayi/A_mag_project/simulations/FPU/'+ str(current_time)  #all files saved in here
    if not os.path.exists(dat_hd):os.mkdir(dat_hd)
    
    #mtf.dump_pickled_data(dat_hd, 'simulation_wave_data', U)
    filename = dat_hd + './test7.hdf5'
    h5py.example_2(filename,U)
    
    
##plotting results###
    #
    ##
    #pm.fft(U,dat_hd)
    #pm.breatherplot(U,dat_hd)
    ##pm.makemovie(U,dat_hd)
    return U

def main():
    N = 10**7
    num_gyro = 16
    simulate(N,num_gyro)
    
    
if __name__ == '__main__':
    main()
