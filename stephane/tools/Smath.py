

from math import *
import numpy as np
import scipy


def cart2pol(x, y):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return(rho, phi)

def pol2cart(rho, phi):
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return(x, y)


def min_local(x,a,w_len=1,window='flat'):
    #find the local minima of x, such that x[i]=np.min(x[i-a:i+a])
    n = len(x)
    indices = []
    valeurs = []
    if w_len>1:
        x=scipy.ndimage.gaussian_filter1d(x,w_len)

    for i,x0 in enumerate(x):
        imin = np.max([0,i-a])
        imax = np.min([i+a,n])

 #       print(x0,np.min(x[imin:imax]))
        
        if x0==np.min(x[imin:imax]):
            #print(imin,imax)
            #print(x0,np.min(x[imin:imax]))
            indices.append(i)
            valeurs.append(x0)
    return indices,valeurs

def max_local(x,a,w_len=1,window='flat'):
    #find the local minima of x, such that x[i]=np.min(x[i-a:i+a])
    n = len(x)
    indices = []
    valeurs = []
    if w_len>1:
        x=scipy.ndimage.gaussian_filter1d(x,w_len)
        #x=np.smooth(x,window_len=w_len,window=window)
        
    for i,x0 in enumerate(x):
        imin = np.max([0,i-a])
        imax = np.min([i+a,n])

        if x0==np.max(x[imin:imax]):
            indices.append(i)
            valeurs.append(x0)
    return indices,valeurs
