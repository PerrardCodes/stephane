# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 15:06:28 2015

@author: stephane
"""

import numpy as np
import stephane.analysis.basics as basics
import stephane.display.graphes as graphes

import stephane.analysis.cdata as cdata
#two assets :
# one point measurement with spatial average along one or two axis
# velocity increments (longitudinal and transverse) for differents spatial distances

def velocity_distribution(M,start,end,display=False):
    #compute the distribution of velocity for Ux, Uy and U for all the individual measurements between start and end
#substract the mean flow in each point

    M = cdata.rm_nan(M,'Ux')
    M = cdata.rm_nan(M,'Uy')
    

    (nx,ny,n)=M.shape()
    nt=end-start

    Ux=np.reshape(M.Ux[:,:,start:end],(nx*ny*nt,))
    Uy=np.reshape(M.Uy[:,:,start:end],(nx*ny*nt,))   
    
    Ux_rms=np.std(Ux)
    Uy_rms=np.std(Uy)
    
    Ux_moy=np.reshape(np.mean(M.Ux[:,:,start:end],axis=2),(nx,ny,1))
    Uy_moy=np.reshape(np.mean(M.Uy[:,:,start:end],axis=2),(nx,ny,1))
    
    Ux_m=np.reshape(np.dot(Ux_moy,np.ones((1,1,nt))),(nx,ny,nt))
    Uy_m=np.reshape(np.dot(Uy_moy,np.ones((1,1,nt))),(nx,ny,nt))
    
#    Ux=np.reshape(M.Ux[:,:,start:end]-Ux_m,(nx*ny*nt,))
#    Uy=np.reshape(M.Uy[:,:,start:end]-Uy_m,(nx*ny*nt,))
    
    Ux=np.reshape(M.Ux[:,:,start:end],(nx*ny*nt,))
    Uy=np.reshape(M.Uy[:,:,start:end],(nx*ny*nt,))
    
    
#    U_s=np.zeros(len(Ux)+len(Uy))
    U_s=np.concatenate((Ux,Uy))
#    U=np.sqrt(Ux**2+Uy**2)
    
    #normalized by the RMS velocity :
    Uxt_rms=np.std(Ux)
    Uyt_rms=np.std(Uy)
    U_rms=np.std(U_s)
    print('RMS velocity : '+str(U_rms)+' m/s')

    mid=(start+end)/2
    
    #Normalisation by the temporal decay function
    Nvec=(M.t[mid]/100)**(-1)
    Nvec = 1
    if display:
        print(max(U_s))
        print(min(U_s))
        
        print(U_s.shape)
        print(Nvec)
      #  graphes.hist(Ux,Nvec,0,100,'o') 
      #  graphes.hist(Uy,Nvec,0,100,'s')
        graphes.hist(U_s,Nvec,fignum=1,num=10**4,label='o') 
         
        title = ''
#        title='Z= '+str(M.param.Zplane)+' mm, t='+str(M.t[mid])+' ms'+', Dt = '+str(nt*M.ft)+' ms'
        graphes.legende('$U_{x,y} (m/s)$','$pdf(U)$',title) 
     #   fields={'Z':'Zplane','t',}
     #   graphes.set_title(M,fields)
    
    return Ux_rms,Uy_rms,Uxt_rms,Uyt_rms

def v_increment(M,start,end,d,p=1,ort='all',fignum=1,normalize=False):
    """
    Compute the distribution of velocity increments, either longitudinal, transverse, or all
    INPUT 
    -----
    M : Mdata object
        with attributes : Ux, Uy
        with method : shape()
    start : int
        start indice
    end : int
        end indice
    d : numpy 1d array
        vector d for computing increments
    p : int
        order of the increments ∂u_p = (u(r+d)^p-u(r)^p)^1/p
    ort : string
        orientation. can be either 'all','trans','long'
    
    """
    #compute the distribution of velocity for Ux, Uy and U for all the individual measurements between start and end    
    (nx,ny,n)=M.shape()
    nt=end-start
    Ux=M.Ux[...,start:end]
    Uy=M.Uy[...,start:end]
    Uz=M.Uz[...,start:end]
    
    
    dim = len(M.shape())
    if dim==3:
        if d[0]>0 and d[1]>0:
            dU_x = (Ux[d[0]:,d[1]:,:]-Ux[:-d[0],:-d[1],:])**p#**(1./p)  #longitudinal component
            dU_y = (Uy[d[0]:,d[1]:,:]-Uy[:-d[0],:-d[1],:])**p#**(1./p)  #transverse component
            dU_y = (Uz[d[0]:,d[1]:,:]-Uz[:-d[0],:-d[1],:])**p#**(1./p)            
        else:
            dU_x = (Ux[d[0]:,...]-Ux[:-d[0],...])**p#**(1./p)
            dU_y = (Uy[d[0]:,...]-Uy[:-d[0],...])**p#**(1./p)   
            dU_z = (Uz[d[0]:,...]-Uz[:-d[0],...])**p#**(1./p)                             
    else:
        print('not implemented')
        
#    U=np.sqrt(Ux**2+Uy**2)

#    graphes.hist(U,1,100,'k^') 
    graphes.hist(dU_x,fignum=fignum,num=10**3,label='ro',log=True) 
    graphes.hist(dU_y,fignum=fignum,num=10**3,label='bs',log=True)
    graphes.hist(dU_z,fignum=fignum,num=10**3,label='m^',log=True)
    
    mid=(start+end)/2
   # title='Z= '+str(M.param.Zplane)+' mm, t='+str(M.t[mid])+' ms'+', Dt = '+str(nt)
    figs = {}
    figs.update(graphes.legende('$dU_{x,y}$','rho(U)','D = '+str(d[0])))
    
    return figs
    
     
def velocity_profile(M,xlines,ylines,display=True,start=0,end=10000,label='k^'):
    nx,ny,nt=M.shape()
    nt=min(nt,end)
    U=np.sqrt(M.Ux[:,:,start:nt]**2+M.Uy[:,:,start:nt]**2)
    label=['k^','rx','bo']
    
    Dt=10
    
    t=M.t[start+Dt:nt-Dt]
    Ut=[]
    for i in ylines:
        for j in xlines:
            Ut.append(basics.smooth(U[i,j],Dt))#[np.mean(S.Uy[i,j,k-Dt:k+Dt]) for k in range(Dt,nt-Dt)]

#            std_U=[np.std(U[i,j,k-Dt:k+Dt]) for k in range(Dt,nt-Dt)]
            if display:
                graphes.graph(t,Ut[-1])
                graphes.legende('t (ms)','V (m/s)','')
    #return a list of time series, for each element in xlines and ylines 
    return t,Ut

def velocity_profile_xy(S,xlines,ylines,display=False,label='k^'):
    nx,ny,nt=S.shape()
    label=['k^','rx','bo']
    
    t=S.t
    Dt=5
    Uxt=[]
    Uyt=[]
    for i in ylines:
        for j in xlines:
          #  std_Ux=[np.std(S.Ux[i,j,k-Dt:k+Dt]) for k in range(Dt,nt-Dt)]
           # std_Uy=[np.std(S.Uy[i,j,k-Dt:k+Dt]) for k in range(Dt,nt-Dt)]

            Uxt.append(basics.smooth(S.Ux[i,j],Dt))#(-1)**i*(-1)**j*    [(-1)**i*(-1)**j*np.mean(S.Ux[i,j,k-Dt:k+Dt]) for k in range(Dt,nt-Dt)]
            Uyt.append(basics.smooth(S.Uy[i,j],Dt))#[np.mean(S.Uy[i,j,k-Dt:k+Dt]) for k in range(Dt,nt-Dt)]

            if display:
#                plt.subplot(211)
                graphes.graph(t[Dt:-Dt],Uxt[-1])#,std_Ux)
                graphes.legende('t (ms)','V (m/s)','Ux')
                
 #               plt.subplot(212)
                graphes.graph(t[Dt:-Dt],Uyt[-1])#,std_Uy)
                graphes.legende('t (ms)','V (m/s)','Uy')
                
    return t,Uxt,Uyt
    