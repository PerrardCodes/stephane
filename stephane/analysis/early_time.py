# -*- coding: utf-8 -*-
"""
Created on Tue Oct 20 10:17:02 2015

@author: stephane
"""

#Focus mainly on the generation process, and the first 5s of the movies (isotropy, homogeneity, vertical spreading)


# isotropy : 
# for each time, compute the distribution of angles (!) in horizontal plane 
# at a very earlier time, should be peaked along x and y directions. Then spread in box directions. Compute that from 5000 or 10000fps movies

import numpy as np
import stephane.display.graphes as graphes

def isotropy(M,label='k^--',display=True,fignum=1):
    step = 1
    tl = M.t[0:None:step]

    N = 50
    display_part = False
    
    Anisotropy = np.zeros(len(tl))
    Meanflow = np.zeros(len(tl))
    
    for i,t in enumerate(tl):
        print(i*100/len(tl))
        rho,Phi = angles(M,i)            
        
        theta,U_moy,U_rms = angular_distribution(M,i)
#        t,U_moy,U_rms = time_window_distribution(M,i,Dt=40)

        if display_part:
            graphes.hist(Phi,fignum=1,num=N)
            graphes.legende('Phi','PDF','')
        
            graphes.graph(theta,U_moy,fignum=3,label='k^')
            graphes.legende('$\theta$','$U^p$','Angular fluctation distribution')

            graphes.graph(theta,U_rms,fignum=4,label='ro')
            graphes.legende('$\theta$','$U^p$','Angular average flow')
        
        Anisotropy[i] = np.std(U_rms)/np.nanmean(U_rms)
        Meanflow[i] = np.std(U_moy)/np.nanmean(U_rms)

    graphes.semilogx(tl,Anisotropy,label='ro',fignum=fignum,subplot=(1,2,1))
    graphes.legende('Time (s)','I','Anisotropy'+graphes.set_title(M))
    graphes.set_axes(10**-2,10**4,0,2)
    graphes.semilogx(tl,Meanflow,label='k^',fignum=fignum,subplot=(1,2,2))
    graphes.legende('Time (s)','<U>','Average flow')
    graphes.set_axes(10**-2,10**4,0,4)
    
    #    input()

#homogeneity makes sense only on statistical average between different recordings.
# the spreading might be easier to look at   
        
def angles(M,i):
    edg = 5    
    Ux = M.Ux[edg:-edg,edg:-edg,i]
    Uy = M.Uy[edg:-edg,edg:-edg,i]  
    
    rho,phi = cart2pol(Ux,Uy)
    return rho,phi
   
def angular_distribution(M,i,p=1,N=100):
    dtheta = np.pi/N
    theta_list = np.arange(0,np.pi,dtheta)
    
    U_moy = np.zeros(N)
    U_rms = np.zeros(N)
    
    for j,theta in enumerate(theta_list):
        U_moy[j],U_rms[j] = U_average_angle(M,i,theta,p=p)
    
    return theta_list,U_moy,U_rms
              
def U_average_angle(M,i,theta,p=1):
    #compute the rms velocity (energy) in the given direction theta
    Ux = M.Ux[:,:,i]
    Uy = M.Uy[:,:,i]
    
    U_moy = np.nanmean(Ux*np.cos(theta) + Uy*np.sin(theta))
    U_rms = np.nanstd(Ux*np.cos(theta) + Uy*np.sin(theta))
    
    return U_moy,U_rms
    
def time_distribution(M,i,Dt=40):
    #at each point, conmpute the average flow and the Urms velocity on a time scale Dt
    # Dt should be greater (!) than the correlation time
    # pointless : can only be used by averaging several movies
    pass

def cart2pol(x, y):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return(rho, phi)

def pol2cart(rho, phi):
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return(x, y)