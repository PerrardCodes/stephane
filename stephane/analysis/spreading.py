# -*- coding: utf-8 -*-
"""
Created on Wed Oct 21 15:34:54 2015

@author: stephane
"""

import stephane.display.graphes as graphes
import numpy as np
import stephane.analysis.basics as basics
import stephane.analysis.vgradient as vgradient

import stephane.pprocess.check_piv as check

import scipy.interpolate as interp

def spreading(Mlist,logscale=False):
    labels = ['k^','ro','bs']
    M = Mlist[0]
    nx,ny,nt = M.shape()
    
    #for i,M in enumerate(Mlist):
    #    vertical_spreading(M,10,Dt=20,label=labels[i])
    print(M.shape())
    frame = range(10,350,10)
    Dt = np.asarray([M.t[frame[i]]-M.t[frame[i-1]] for i in range(1,len(frame))])
   # print(Dt)
#    Ux_adv = np.zeros(nx)
#    Uy_adv = np.zeros(ny)    

    D = 5*10**2
    z_init,U_init = initial_shape(offset=140)
    
    t_init = 0.026

    for i,t0 in enumerate(frame[:-1]):
        z,Ux_moy,Uy_moy,Ux_rms,Uy_rms = profile_average(Mlist,t0,Dt=10,display=False)
        
        print(M.t[t0])
        z_dif,U_dif = solve_diffusion(-z_init,U_init,D,M.t[t0]-t_init)
        
       # print(U_dif)
        
        if i==0:
            fx = interp.splrep(z, Ux_rms,s=2)# kind='cubic')
            fy = interp.splrep(z, Uy_rms,s=2)# interp1d(z, Ux_rms, kind='cubic')
                
            Ux_adv = interp.splev(z, fx)
            Uy_adv = interp.splev(z, fy)

        figs = display_profile(Mlist[0],i,t0,z,Ux_moy,Uy_moy,Ux_rms,Uy_rms,Dt=Dt[i],fig=i,logscale=logscale)
   #     Ux_adv,Uy_adv = advection(Dt[i],i,z,Ux_moy,Uy_moy,Ux_adv,Uy_adv,display=True)
              
        graphes.semilogx(U_dif,z_dif,fignum=i,label='r-')
        graphes.set_axis(10**0,10**4,-300,-120)        
        
        
        graphes.save_figs(figs,savedir='./Spreading/Diffusion_noNoise/',suffix='',prefix='2015_12_28_front_'+str(t0),frmt='png')
        
    figs = {}    
    graphes.set_fig(1)
    figs.update(graphes.legende('$t$ (s)','$Front position z$ (mm)',''))
    graphes.set_fig(2)    
    figs.update(graphes.legende('$t$ (s)','$Front position z (log)$ (mm)',''))
    graphes.set_fig(3)    
    figs.update(graphes.legende('$t$ (s)','$U_{rms}$ (mm/s)',''))

    graphes.save_figs(figs,savedir='./Spreading/',suffix='',prefix='2015_12_28_all_')    
#    graphes.save_figs(figs,savedir='./Spreading/',suffix='',prefix='2015_12_28_'+str(i)+'_')

    
def advection(Dt,i,z,Ux_moy,Uy_moy,Ux_rms,Uy_rms,display=True):
    """
    Advect the current profile by the instantaneous x-averaged Uz profile
    """
    #first spline the vertical profile
    #then,
    fx = interp.splrep(z, Ux_rms,s=0)# kind='cubic')
    fy = interp.splrep(z, Uy_rms,s=0)# interp1d(z, Ux_rms, kind='cubic')
    
    dz = np.diff(z[0:2])[0]
    dzi = dz/10.
    zi = np.arange(min(z)-dzi,max(z)+dzi,dzi)
    
    dU_x = interp.splev(z, fx, der=1)
    dU_y = interp.splev(z, fy, der=1)
    
    V = Ux_moy
   # dt = M.t[t0+1]-M.t[t0]
    
    print(np.mean(Ux_rms))
    print(np.max(Dt * V * dU_x))
    print(np.min(Dt * V * dU_x))
        
    Ux_adv = Ux_rms - Dt * V * dU_x
    Uy_adv = Uy_rms - Dt * V * dU_y
    
    if display:
        graphes.semilogx(Ux_adv,z,fignum=i,label='r-')
        graphes.semilogx(Uy_adv,z,fignum=i,label='r-')
    
    return Ux_adv,Uy_adv
    
def vertical_spreading(M,N,Dt=10,label='k^'):
    """
    Compute the vertical profile of RMS velocity as a function of time
    
    INPUT
    -----
    
    OUTPUT
    -----
    figs : dict
        key correspond to figure index, value to their standarzied filename
    """
    n=160
    ti=50
    
    figs={}
    
    z0 = - M.param.stroke #!!!
    z_1=np.zeros(n)
    z_2=np.zeros(n)
    
    t=np.zeros(n)
    E=np.zeros(n)
     
    indices = [(i+1)*N+ti for i in range(n)]
    for i,t0 in enumerate(indices):
        z_1[i],z_2[i],E[i],std_E = velocity_profile(M,i,t0,Dt=Dt,display=True)         
        t[i] = M.t[t0]
    #average on horizontal line (if rotated, correspond to vertical line)
#compute the vertical RMS velocity profile

    fig=1
    graphes.set_fig(fig)
    graphes.graph(t,z_1-z0,fignum=fig,label=label)
    graphes.graph(t,z_2-z0,fignum=fig,label=label)
    figs.update(graphes.legende('$t$ (s)','$Front position z$ (mm)',''))

    fig=2
    graphes.set_fig(fig)
    graphes.graphloglog(t,np.abs(z0-z_1),fignum=fig,label=label)
    graphes.graphloglog(t,np.abs(z0-z_2),fignum=fig,label=label)
    
    graphes.graphloglog(t,np.power(t,0.5)*10**2,fignum=fig,label='r--')
    figs.update(graphes.legende('$t$ (s)','$Front position z (log)$ (mm)',''))

    fig=3
    graphes.set_fig(fig)
    graphes.graphloglog(t,E,fignum=fig,label=label)
    figs.update(graphes.legende('$t$ (s)','$U_{rms}$ (mm/s)',''))
    
    t_min = 6*10**-2
    t_max = 4*10**-1
    indices = np.where(np.logical_and(t<t_max,t>t_min))
    t_part = t[indices]
    z_part = np.abs(z0-z_1[indices])

    P = np.polyfit(np.log(t_part/t_min),np.log(z_part),1)
    C = np.exp(P[1])
    nu = C**2 / t_min
    print(P[0],C)
    print('Effective diffusion coefficient : '+str(nu)+' mm^2/s')
    
    graphes.graphloglog(t,C*np.power(t/t_min,P[0]),fignum=2,label='r-')
    figs.update(graphes.legende('$t$ (s)','$Front position z (log)$ (mm)',''))

    graphes.save_figs(figs,savedir='./Spreading/',suffix='',prefix='2015_12_28_front_both_C')    

    return figs
    
            
def profile_average(Mlist,t0,Dt=10,display=False):
    nx,ny,nt = Mlist[0].shape() #assume for now that all the experiments have exactly the same data dimensions
    Ux = np.zeros((ny,len(Mlist)))
    Uy = np.zeros((ny,len(Mlist)))
    Ux_r = np.zeros((ny,len(Mlist)))
    Uy_r = np.zeros((ny,len(Mlist)))

    for i,M in enumerate(Mlist):
        z,Ux[:,i],Uy[:,i]=basics.profile_average(M,t0,t0+Dt,p=1)
        z,Ux_r[:,i],Uy_r[:,i]=basics.profile_average(M,t0,t0+Dt,p=2)
    
    Ux_moy = np.nanmean(Ux,axis=1)
    Uy_moy = np.nanmean(Uy,axis=1)
    Ux_rms = np.nanmean(Ux_r,axis=1)
    Uy_rms = np.nanmean(Uy_r,axis=1)
    
    return z,Ux_moy,Uy_moy,Ux_rms,Uy_rms
    
def fit_profile(M,t0,z,U,p=9,log=True):
    
    U_noise_low,U_noise_high = check.bounds(M,t0)    
    
    if log:
        Y = np.log10(U)
    else:
        Y = U
        
    P = np.polyfit(z,Y,p)#-U_noise_low,p)
    U_th = np.polyval(P,z)
    
    if log:
        U_th=np.power(U_th,10)
        
    return z,U_th
        
def display_profile(M,i,t0,z,Ux_moy,Uy_moy,Ux_rms,Uy_rms,Dt=1,fig=1,logscale=True):
    t = M.t[t0]
    U_noise_low,U_noise_high = check.bounds(M,t0)    
    graphes.set_fig(fig)
    title = 't='+str(int(t*1000))+' ms'+'Urms_zprofile'

    Dir=M.fileDir+'Velocity_distribution_log_M_2015_12_28_Meanfield'+'/'
    
    z,Ux_th = fit_profile(M,t0,z,Ux_rms,p=9,log=True)
    z,Uy_th = fit_profile(M,t0,z,Uy_rms,p=9,log=True)
    
  #  Ux_adv,Uy_adv = advection(Dt,z,Ux_moy,Uy_moy,Ux_rms,Uy_rms)
  #  print(Ux_adv/Ux_rms)
#    print(Ux_adv)
    
    
    figs={}
    if logscale:
     #   graphes.semilogx(Ux_rms,z,fignum=0,label='bo--')
        graphes.semilogx(Uy_rms,z,fignum=i,label='k^--')
                
        graphes.semilogx([U_noise_low,U_noise_low],[-400,-100],fignum=i,label='r--')
        graphes.semilogx([U_noise_high,U_noise_high],[-400,-100],fignum=i,label='r--')

     #   graphes.semilogx(np.sqrt(np.power(Ux_moy,2)),z,fignum=i,label='b+--')
     #   graphes.semilogx(np.sqrt(np.power(Uy_moy,2)),z,fignum=i,label='c ^--')
        
        graphes.set_axis(10**0,10**4,-300,-120)        
        figs.update(graphes.legende('$U_{rms} (t/t_0)$','$z$ (m)',''))
        file=graphes.set_title(M,title)
        filename=Dir+file
        
    else:
        graphes.graph(Ux_rms,z,fignum=0,label='bo--')
        graphes.graph(Uy_rms,z,fignum=i,label='k^--')
        
     #   graphes.semilogx(Ux_th,z,fignum=i,label='r-')
     #   graphes.semilogx(Uy_th,z,fignum=i,label='r-')
        graphes.graph([U_noise_low,U_noise_low],[-400,-100],fignum=i,label='r--')
        graphes.graph([U_noise_high,U_noise_high],[-400,-100],fignum=i,label='r--')

        graphes.graph(np.sqrt(np.power(Ux_moy,2)),z,fignum=i,label='b+--')
        graphes.graph(np.sqrt(np.power(Uy_moy,2)),z,fignum=i,label='c ^--')
        
        graphes.set_axis(0,2.5*10**3,-300,-120)        
        figs.update(graphes.legende('$U_{rms} (t/t_0)$','$z$ (m)',''))
        file=graphes.set_title(M,title)
        filename=Dir+file
        
        graphes.save_figs(figs,savedir='./Spreading/Stat_average_lin/',suffix='',prefix='2015_12_28_front_'+str(t0),frmt='png')
            
    return figs


def front():
    pass
    
def front_time(M,i,t0,Dt=10,display=True):
    """
    Compute the time of propagation of a front.
    From the temporal signal of U_rms in one horizontal
    INPUT
    -----
    
    OUTPUT
    -----
    """
    pass
    
def velocity_profile(M,t0,Dt=10,display=False):
    """
    Comnpute the velocity profile of U_rms along the vertical direction.
    both components of U_rms are obtained by an average over a time Dt and a horizontal line in space
    
    INPUT
    -----
    M : Mdata object
        data
    i : int 
        index in time
    t0 : 
    OUTPUT
    -----
    """
    mid=t0+Dt/2
 #   t[i]=M.t[mid]   
#    rect = [16,len(M.x[])]
#    basics.crop(M,rect)
    z,Ux_moy,Uy_moy=basics.profile_average(M,t0,t0+Dt,p=1)
    z,Ux_rms,Uy_rms=basics.profile_average(M,t0,t0+Dt,p=2)
    
    z_min = -190
    z_max = -170 
    indices = np.where(np.logical_and(z<z_max,z>z_min))
    E0,std_E0 = (np.median(Uy_rms[indices]),np.std(Uy_rms[indices]))
    
    cut = 16
    z_min = -300
    z_max = -150
    indices = np.where(np.logical_and(z<z_max,z>z_min))
        
    Ux_rms = Ux_rms[indices]
    Uy_rms = Uy_rms[indices]
    Ux_moy = Ux_moy[indices]
    Uy_moy = Uy_moy[indices]
    z = z[indices]

    z_Ux,z_Uy,Ux_norm,Uy_norm = half_width(z,Ux_rms,Uy_rms,E0)
    
    return z_Ux,z_Uy,E0,std_E0

    
def half_width(z,Ux_rms,Uy_rms,E0):
    
    z0=-170
    indice = np.argmin(np.abs(z-z0))
    #    print(z,Ux_rms,Uy_rms)
    Ux_norm=Ux_rms/Ux_rms[indice]#*M.t[mid]/0.1 
    Uy_norm=Uy_rms/Uy_rms[indice]#*M.t[mid]/0.1
#    print('Ux_0 : '+str(Ux_rms[indice]))
    
    thres = 0.5
    ind = np.argmin(np.abs(Ux_norm-thres))
    z_Ux = z[ind]
    
    ind = np.argmin(np.abs(Uy_norm-thres))
    z_Uy = z[ind]
#    print('z_Ux : '+str(z_Ux))
#    print('z_Uy : '+str(z_Uy))
    
    return z_Ux,z_Uy,Ux_norm,Uy_norm

def signal_noise(M,z,mid):

    U_noise_low,U_noise_high = check.bounds(M,mid)

    #find the position where the fluctuating velocity is 4time greater than the noise level
    indice = np.argmin(np.abs(4*U_noise_low-Uy_rms))
    z_front = z[indice]

    
    #catch 3 points before and after the cut of
    step = 5
    imin = np.max([0,indice-step])
    imax = np.min([len(z)-1,indice+step])
    
    z_part = z[imin:imax]
    Ux_part = Ux_rms[imin:imax]
    Uy_part = Uy_rms[imin:imax]
    
#    interpolate by a linear function to find the exact position of the front :    
    Uy_part-4*U_noise_low
    P = np.polyfit(z_part-z_front,Uy_part-4*U_noise_low,1)
    
    if np.abs(P[1])<np.abs(np.mean(np.diff(z_part))):
        z_front += P[1]
    else:
        pass
      
def front(M):
    """
    From a temporal measurement of a turbulent spreading situation,
    compute the position of the front.
    The algorithm is based on the following steps :
    Start with an horionzal line on the top.
    At each time step, evaluates if the front has been moving downward.
    Detection based on ? threshold of detection ? energy amplitude ?? 
    Update the position of the front
    """
    
    X,Y = graphes.get_axis_coord(M)
    nx,ny,nt = M.shape()
    
    Z_front = np.zeros(nx)+75  #the front is propagating downward
    
    #Z_front = np.max(Y,axis=1) #?? should correspond to the X axis ???
    
    frame = 100
    #energy profile along the front :
    
#    print(len(Z_front))
    field = 'E'
    if not hasattr(M,field):
        M,field = vgradient.compute(M,field,filter=True)
    
    test_front(M,Z_front,frame)
    
def test_front(M,Z,frame):
    nx,ny,nt = M.shape()
    Z_test = Z + np.ones(nx) 
            
    E_front = np.asarray([M.E[i,j,frame] for i,j in enumerate(Z)])
    Emoy = np.nanmean(E_front)
    print(Emoy)
    print(np.nanstd(E_front))
    
    thres = 0.5
    for i,j in enumerate(Z_test):   
        if M.E[i,j,frame]>Emoy*thres:
            Z[i]+=1
        
    print(Z_test-Z+1)
    

def diffusion(M):
    
    z,U = initial_shape(offset=0)
    
    dt = 0.05
#    print(M.t)
    times = M.t[100:1000:100] #np.arange(dt,0.5,dt)
    
    D = 10**4
    for t in times:
        solve_diffusion(z,U,D,t)
    
    z,U = initial_shape(offset=0)
    
    
def initial_shape(offset=0):
    zmin = -500+offset
    zmax = 500+offset
    
#    Dz = 45
    
    dz = 1
    z = np.arange(zmin,zmax,dz)

    Dz = 200    
    U0 = 2.5*10**3
    delta = 10.
    
    U1 = U0*np.exp(-Dz/delta)
    #delta = Dz / np.log(U0/U1)

    z0 = 0+offset
    z1 = Dz+offset
    
    U = np.zeros(len(z))
    for i,zi in enumerate(z):
        if zi<z0:
            U[i] = U0
        if zi>=z0 and zi<z1:
            U[i] = U0*np.exp((z0-zi)/delta)
        if zi>=z1:
            U[i] = U0*np.exp(-Dz/delta)

    graphes.semilogy(z,U,label='r')
    graphes.legende('z (mm)','U (mm/s)','')
    graphes.set_axis(-50,150,10**0,5*10**3)
    
    return z,U
    #front format : f 1d np array gives the vertical position for each X position.
    
def solve_diffusion(z,U,D,t):
    
    U_th = 1/np.sqrt(4*np.pi*D*t)*np.exp(-np.power(z,2)/(4*D*t))

    U_sol = np.zeros(len(U))
    
    for i,zi in enumerate(z):
        U_th = 1/np.sqrt(4*np.pi*D*t)*np.exp(-np.power(zi-z,2)/(4*D*t))        
        U_sol[i] = np.sum(U_th*U)
        
 #   print(U_sol)
   # graphes.semilogy(z,U_sol,fignum=1,label='k')
   # graphes.set_axis(-50,150,10**0,5*10**3)
    
    return z,U_sol
     #    a=-2.4
     #    b=0.12        
       #      U_spatial=np.concatenate((a*z[z<0]+b,0*z[z>0]+b))
                        
         #find the z position where Ux_rms ~ 0.2 m/s
             # plot this z position as a function of time
        
     #    thres=0.3 
     #    indices=np.where(Ux_norm<thres)
     #    if len(indices)>1:
     #            print(indices)
     #            print(len(indices))
      #       indice=indices[0][0]
      #   else:
       #      indice=np.argmin(np.abs(Ux_norm-thres)) #Uy correspond to the Ux component for angle=90    
          #   print(i,indice)
     
      
      