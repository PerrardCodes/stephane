# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 11:50:34 2015

@author: stephane
"""

import numpy as np
import stephane.display.graphes as graphes
import stephane.analysis.cdata as cdata
import stephane.manager.access as access
#import time_step_sample

def decay(M,field,M_moy=None,display=False,label='',fignum=1,compute=True,fluctuations=False,log=True):
    """
    Compute the spatial-averaged time evolution of kinetic energy of a flow field
    
    INPUT
        M : Mdata object
            Data set containing a 2D dimensionnal flow (files Ux and Uy) as a function of time 
        display : boolean. defaulf value False
            To display the result of the computation
        label : standard matplotlib label format
        fignum : standard graphes.py fignum format
            (-1 : current figure, 0: clear figure, N(>1): set figure N)
    OUTPUT
        tf : 1d numpy array
            time axis
        Uf : 1d numpy array.
            Spatial average kinetic energy
    """
    Y = M.get(field)
    
    # compute decay of turbulent kinetic energy from the total energy.
    # to compute from a fit of the energy cascade go to Fourier.py
    if M_moy is not None:
        Ux = M.Ux - M_moy.Ux
        Uy = M.Uy - M_moy.Uy
    else:
        Ux = M.Ux
        Uy = M.Uy
        
    if compute and field=='E':
        Y = (Ux**2+Uy**2)
    else:
        Y = access.get_all(M,field)
        
#    if fluctuations:
#        Y = access.get_all(M,field) - (M.Ux**2+M.Uy**2)
    
    nx,ny,nt=Y.shape
    Y_t=np.reshape(Y,(nx*ny,nt))    
    #retake the time instants from the cinefile itself
#    times=time_step_sample.get_cine_time(M.Sdata.fileCine,False)    
    #current : directly take the attribut t of M :
    times = M.t
    
    t=[]
    U_rms=[]
    #time of the free fall
    t0=0.#.53#0;550/M.ft
    for i in range(nt):
#        print(i)
#        print(Dt[i],np.nanmedian(E_t[:,i]/Dt[i]**2,axis=0))
        
   #     U_rms.append(np.nanmedian(Y_t[:,i],axis=0))
        
        U_rms.append(np.sqrt(np.nanmedian(np.power(Y_t[:,i],2),axis=0)))
        
#        U_rms.append(np.nanmean(E_t[:,i]/Dt[i]**2,axis=0))
      #  t.append(times[int(M.im_index[i])]-t0)
        t.append(times[i]-t0)
        
    tf=np.asarray(t)
    Uf=np.asarray(U_rms)
    
    figs = {}
    if display:
#        figs = {}
    #    graphes.graph(tf,Uf,fignum=fignum,label=label)
        if log:
            graphes.graphloglog(tf,Uf,fignum=fignum,label=label)
            graphes.set_axis(10**-2,10**2,10**-1,10**6)
 #       graphes.graphloglog([10**-2,10**2],[10**-1,10**4],fignum=-2,label='r--')
     #  graphes.graphloglog([10**-1,10**3],[5*10**7,5*10**-1],fignum=-2,label='r--')
        else:
           graphes.graph(tf,Uf,fignum=fignum,label=label)
            
        figs.update(graphes.legende('$t$ (s)','$<E>_{x,y}$ (mm^2/s^2)',graphes.title(M)))
        return figs
    else: 
        return tf,Uf
    
def correlation_length():
    #compute the correlation length
    #for every point, look at any point at a given distance d. Compute sum(U(x)*U(x+r cos theta))
    # -> now see corr.py    
    u[i]*(U[j]-Umoy)
    
    
def mean_flow(M_log):
    #compute the mean flow by averaging the flow on the 5 movies and 10 times steps (more or less depending on the time ratio !!)
    Ux=[]
    Uy=[]

    M_ref = M_log[0]   
#    X=M_ref.x
#    Y=M_ref.y

    t=M_ref.t
    Etot = []    
    for tp in range(len(t)):
        Umoy,Vmoy=average(M_log,tp)
        
        E = np.sqrt(Umoy**2+Vmoy**2)
        
        Etot.append(np.nanmean(E))
        print(np.nanmean(E))
#        graphes.color_plot(X,Y,E,fignum=3)

        Ux.append(Umoy)
        Uy.append(Umoy)
 
    graphes.graph(t,Etot,fignum=-1,label='k^')
    graphes.legende('$t$ (s)','$\sqrt{<E_c>}$ (mm/s)','')
       
    graphes.graphloglog(t,Etot,fignum=-1,label='k^')
    graphes.legende('$t$ (s)','$\sqrt{<E_c>}$ (mm/s)','')
       
    return Ux,Uy
    
def fluctuation():
    pass
    #same as mean flow, but for the standard deviation between one movie to another.
    #compute a typical length scale associated to the flow fluctuation. Correlation function between two points ? (as a function of distance)

def main():  
#    mean_flow(M_new+M_log)
    labels = ['k^','ro','bp','c8','g*']

    for i,label in enumerate(labels):
        M = Mlist2[i]
        t_d,E = decay(M,label=label,display=True)
        input()
        
        #2015_09_22 : good indices, 0,1,3 (?). bad one : 2,4
#    for M in M_zoom_zoom:
#        decay(M)
#    for M in M_zoom:
#        decay(M)
        
#main()