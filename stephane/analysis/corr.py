# -*- coding: utf-8 -*-
"""
Created on Thu Jul 30 14:35:11 2015

@author: stephane
"""


import numpy as np
import math
import os.path
import scipy.optimize
import stat

import stephane.tools.rw_data as rw_data
import stephane.tools.Smath as Smath
import stephane.tools.fitting as fitting


import stephane.display.graphes as graphes
import stephane.analysis.cdata as cdata

#import stephane.analysis.decay as decay
import stephane.analysis.Fourier as Fourier
import stephane.analysis.statP as statP
#import fit

import stephane.manager.access as access
    
def time_correlation(Mlist,indices=None,display=False):
    """
    Compute the spatial averaged time of velocity autocorrelation
    Velocity autocorrelation functions in time are fitted by an exponential in time.
    Typical time tc gives the time correlation
    INPUT
    -----
    Mlist : list of Mdata
    indices : list of int
        indices of Mlist elements to process. default value process all the elements
    display : bool
        default value False
    OUTPUT
    -----
    """
    
    if indices is None:
        indices = range(len(Mlist))
     
    labels = ['k^','ro','bp','c8','g*']
    
    for i,indice in enumerate(indices):
        label = labels[i]
        
        M = Mlist[indice]
        tf,tau=compute_Ct(M,display=False,label='ko',fignum=1)

        graphes.graphloglog(tf,tau,fignum=9,label=label)
        graphes.legende('$t (s)$','$\tau (s)$','')
    
        #compute from the 
       # t_d,E = decay.decay(M,label=label)
        t_d,E = Fourier.display_fft_vs_t(M,'1d',Dt=50,label=label)
        
        Ef = np.zeros(len(tf))
        for i,t in enumerate(tf):
            j=np.argmin(abs(t_d-t))
          #  print(str(j)+ ' : '+str(E[j]) + ", " + str(tau[i]))
            Ef[i] = E[j]
            
        
        graphes.graphloglog(Ef,tau,fignum=10,label=label)
        graphes.legende('$E (m^2/s^2)$','$\tau (s)$','')
        

def spatial_correlation(M,compute=True,rootdir='Corr_functions',label='k^',fignum=1,display=True,save=False):
    """
    Compute the spatial correlation function 
        or display the correlation length as a function of time 
        save the correlations function in txt files
    INPUT
    -----
    M : Mdata object
    compute : bool
        default value : True
        if True, the correlations functions are computed, and save in txt files 
        if False, load the previously computed correlations functions and display the correlation length as a function of timme
    rootdir : string
        subdirectory name for saving the corr functions. Base directory correspond to the location of the associated dataset.
    OUTPUT
    -----
    None
    """
    
    if compute:
        # chose randomly the pair of indices that will be used for computing the spatial correlation function
        dlist=range(int(max(M.shape())/2))
        indices={}
        N = 10**3
        for i,d in enumerate(dlist):
            indices[i]=d_2pts_rand(M.Ux[:,:,0],d,N)   
        print('Pair of indices computed')
        (ny,nx,nt) = M.shape()
            
      #  tref,d,Cxx,Cyy,Cxy,CEE = correlation_functions(M,dlist,indices,Dt=Dt)
      #  print('Correlation functions computed')
        
        step = 100
        Dt = 20
        tref,d,Cxx,Cyy,Cxy,CEE = display_corr_vs_t(M,dlist,indices,step=step, Dt=Dt,label='-',display=True)
    
        #How to save matrices in a single txt file ?    
    #    filename = os.path.dirname(M.filename) + '/Corr_functions/' + name + '.txt'
        axes = ['xx','yy','xy','CEE']
        keys = ['d','t'] + ['Corr_'+p for p in axes]
    
        Corr_functions = [Cxx,Cyy,Cxy,CEE]    
        List_info = [d,tref] + Corr_functions
        
        if save:
            name = 'Corr_spatial_'+M.Id.get_id()
            filename = os.path.dirname(M.filename) + '/Corr_functions/' + name + '.txt'
            print(filename)
            rw_data.write_matrix(filename,keys,List_info)
       # for C in [Cxx,Cyy,Cxy,CEE]:
       #     rw_data.write_matrix(tref,D,C)
    
            print('Correlation functions saved in '+filename)
    
    else:
        #try first to load the data rfrom the Corr_functions directory
        Dir = os.path.dirname(M.filename) + '/' + rootdir 
        filename = Dir + '/Corr_spatial_' + M.Id.get_id() + '.txt'  #print(filename)
        Header,Data,axes = rw_data.read_matrix(filename,Hdelimiter='\t',Ddelimiter='\t')          
        print(Data.keys())
        nd = axes['d']
        nt = axes['t']
        
        print((nd,nt))    
        Dt = 20
        for key in Data.keys():
            Data[key]=np.reshape(np.asarray(Data[key]),(nd,nt))  
#            Data[key]=cdata.smooth(Data[key],Dt)
        dimensions = M.shape()
        tlist = range(Dt,dimensions[2]-3*Dt,1)
        lc = np.zeros(len(tlist))
        tf = np.zeros(len(tlist))
        
        display_part = True
        for i,t in enumerate(tlist):
            d=Data['d'][:,i]
            
            key = 'Corr_xx'
            Cd=Data[key][:,i]/Data[key][0,i]
            
            popt = fitting.fit(fitting.exp,d[0:10],Cd[0:10])
                    
            tf[i] = M.t[t]
            lc[i]=-1/popt[0]
         #   print(str(M.t[t]) + ' : ' +str(lc[i]))
            if display_part:
                if i%100==0:
                    graphes.set_fig(1)
                    graphes.graph(d,Cd,fignum=0,label='ko')
                    graphes.graph(d,exp(d,-1/lc[i]),label='r')
                 #   graphes.graph(d,parabola(d,-1/lc[i])+1,label='r')
                    graphes.legende('$d (mm)$','$C_d$','')
                    graphes.set_axes(0,3,-1,1.1)
                    #input()
        
        cdata.rm_nans([lc],d=2)
                
        if display:    
            graphes.graphloglog(tf,lc,fignum=fignum,label=label)
            graphes.legende('$t (s)$','$d (mm)$',graphes.title(M))

def correlation_functions(M,dlist,indices,Dt=1):
    """
    compute two points correlation functions
        various axis (xx, yy, xy, E)
        ?? add longitudinal/transverse correlation functions ?
        compute corr function from a list of pair indices previously generated, corresponding to a list of distances
    INPUT
    -----
        M : Mdata object to be processed
        dlist : int array. Each element correspond to a given distance
        indices : list of tuple containing pair of indices 
        Dt : int, default value : 1
            time window for averaging over time
    OUTPUT
    -----
        tref : np array
            time axis
        d : 1d np array
            distance axis
        Cxx : 2d np array
            autocorrelation function along xx
        Cyy : 2d np array
            autocorrelation function along xx
        Cxy : 2d np array
            crosscorrelation function between x and y
        CEE : 2d np array
            autocorrelation function of energy E
    """    
   # tref = [M.t[1000:1010]]
    tref = M.t[50:60]
    
    nt = len(tref)#/10
    nd = len(dlist)

    Cxx = np.zeros((nd,nt))
    Cyy = np.zeros((nd,nt))    
    Cxy = np.zeros((nd,nt))
    CEE  = np.zeros((nd,nt))

    for i,t in enumerate(tref):
        print(str(100.*i/len(tref)) + " % ")
      #  print('Time elapsed : '+str(t)+ ' s')
        for d in dlist:
            Cxx[d,i]= corr_v(M,t,indices[d],axes=['Ux','Ux'],p=1)[0]
            Cyy[d,i]= corr_v(M,t,indices[d],axes=['Uy','Uy'],p=1)[0]
            Cxy[d,i]= corr_v(M,t,indices[d],axes=['Ux','Uy'],p=1)[0]
            CEE[d,i] = corr_v(M,t,indices[d],axes=['E','E'],p=1)[0]
            
    dx = M.x[1,0]-M.x[0,0]
    if dx==0:
        dx = M.x[0,1]-M.x[0,0]
    fx = M.fx/dx
#        fx=M.fx/max([M.x[0,1]-M.x[0,0],)
    d = np.asarray(dlist)*fx
 
 #Smoothing functions           
#    Cxx = cdata.smooth(Cxx,Dt)
#    Cyy = cdata.smooth(Cyy,Dt)
#    Cxy = cdata.smooth(Cxy,Dt)
#    CEE = cdata.smooth(CEE,Dt)
    
    if Dt>1:
        tref = tref[0:1-Dt]
    
    return tref,d,Cxx,Cyy,Cxy,CEE
    
def display_corr_vs_t(M,dlist,indices,step=100, Dt=1,label='-',display=False,fignum=1):

    tref,d,Cxx,Cyy,Cxy,CEE = correlation_functions(M,dlist,indices,Dt=Dt)
    
    #Display successive correlations functions
    times = range(0,len(tref)-3*Dt,step)
    times = range(0,len(tref),step)
    
    if display:
        for t in times:
            graphes.graph(d,Cxx[:,t]/Cxx[0,t],fignum=fignum)
            graphes.set_axis(0,max(d),-1,1.5)
            graphes.legende('d (mm)','C_{xx}','')
            
            graphes.graph(d,Cyy[:,t]/Cyy[0,t],fignum=fignum+1)
            graphes.set_axis(0,max(d),-1,1.5)
            graphes.legende('d (mm)','C_{yy}','')
    
            graphes.graph(d,CEE[:,t]/CEE[0,t],fignum=fignum+2)
            graphes.set_axis(0,max(d),-1,1.5)
            graphes.legende('d (m)','C_{E}','')
    
    return tref,d,Cxx,Cyy,Cxy,CEE

    
def correlation_length(d,X):
#from a correlation function, return the correlation length
#try to fit by different models of corr functions ? (exp ?)
    dc=np.asarray(d)    
    Xc=np.asarray(X)
    #method 1
#    indice = np.argmin(np.abs(Xc/Xc[0]-1/math.exp(1)))
#    dcorr = d[indice]
    #method 2
    N=8
    indices=np.arange(N)
    P=np.polyfit(dc[indices],Xc[indices]/Xc[0],1)
    dcorr=-1/P[0]
    
    return dcorr
    
def d_2pts(u,d,epsilon=1):
    #find in u the list of tuple indexes that are separated by a distance d+- epsilon
    #u is a two dimensionnal matrix
    ny,nx=u.shape
    indices={}
#    i0 = ny//2
#    j0 = nx//2
    #compute all the indices of the neighboors.
    # For any other point, just translate the neighboors, and remove everything that are not inside the squared matrix 
    # OR pick up randomly a fixed number of pairs (1000 ?) for each distance    
    # needs a refined way to compute pair of indices (way too long with the current technic !)
    for i in range(ny//2-15,ny//2+15):
    #    print(i*100//ny)
        for j in range(nx//2-15,nx//2+15):
            tata,ind=cdata.neighboors(u,(i,j),b=d+epsilon/2,bmin=d-epsilon/2)
        #    print(ind)
            indices[(i,j)]=ind
#    print(indices)
    return indices
        
def corr_v(M,t,indices,avg=1,axes=['Ux','Ux'],p=1,average=False):
        C=[]
        C_norm=[]
        
        X,Y=chose_axe(M,t,axes)
        
        if average:
            Xmoy,Xstd=statP.average(X)
            Ymoy,Ystd=statP.average(Y)
        else:
            Xmoy = 0
            Xstd = 0
            Ymoy = 0
            Ystd = 0

        for i,j in indices.keys():  
           # print(indices[i,j])
            k,l= indices[i,j]
            Sp = (X[i,j]-Xmoy)**p*(Y[k,l]-Ymoy)**p   #remove the average in space ?   -> remove by default
            C.append(Sp)# for k,l in indices[(i,j)]])  
            
            Sp_norm = (X[i,j]-Xmoy)**(2*p) 
            C_norm.append(Sp_norm)
            # substract the mean flow ? it shouldn't change the result so much, as there is no strong mean flow
            #-> to be checked
            # how to compute the mean flow : at which scale ? box size ? local average ?
            # Cmoy,Cstd=average(C)      
        Cf = statP.average(C)[0]/statP.average(C_norm)[0]
        return Cf
        
def corr_v_tn(M,t,indices,avg=1,p=1):
    """
    Compute the correlation coefficient at time t for a collection of pair indices
    
    INPUT
        M : Mdata object
        t : time indice
        indices : dictionnary containing pair of indices
    """
    C = []
    x = M.x
    y = M.y
    Ux = M.Ux[...,t]
    Uy = M.Uy[...,t]
    
    Xmoy,Xstd=statP.average(Ux)
    Ymoy,Ystd=statP.average(Uy)
                    
    for i,j in indices.keys():  
        k,l = indices[i,j]
        
        #tangent vector
        u = (x[k,l]-x[i,j],y[k,l]-y[i,j])
        theta,r = Smath.cart2pol(u)
        
        Sp = (X[i,j]-Xmoy)**p*(Y[k,l]-Ymoy)**p   #remove the average in space ?
                                
#                for k,l in zip(indices[(i,j)][0],indices[(i,j)][1]):
        C.append(Sp)# for k,l in indices[(i,j)]])    
                    # substract the mean flow ? it shouldn't change the result so much 
                    # but it does change !
                    # at which scale ? box size ? local average ?
                    # 
    #        Cmoy,Cstd=average(C)       
    return statP.average(C)

        
def corr_v_t(Mlist,t,axes=['Ux','Ux'],N=100,p=1,display=False,save=False,label='^',fignum=1):
    """
    compute the correlation function in time around a given time. Average over space (ie along d-1 dimensions of Mlist[axes[0]])
    INPUT
    ------
    M : Mdata object
        must have attributes : t, fields ('Ux','Uy', ...)
        must have method shape()
    t : int
        time index
    axes : 2 elements string list
        attributes of M to be used for the correlation function
    N : int
        Number of frames before and after time t.
    p : int. default value is 1
        power of the fields C_p(Dt) = U1**p * U2**p / <U1**2p >
    display : bool. default value is false
    """
    #    Ct=[]    
    #how to average ?? -> need a mean on several realization ? or just substract the average value (in space) ?
    M = Mlist[0]
    tlist=M.t
   # print("Compute mean values")
    tmin = max(0,t-N)
    tmax = min(len(tlist)-1,t+N)
    
    #print(tmin)
    #print(tmax)
    t_c = np.asarray(M.t)[np.arange(tmin,tmax)]#-M.t[t]
#    print(t_c)
   # print(np.mean(np.diff(t_c)))
    #max number of time step 
  #  print("Compute mean values 2")  
    #Compute the average flow
    Xm = []
    Ym = []
    for M in Mlist:
        Xref,Yref=chose_axe(M,t,axes)
        Xm.append(Xref)
        Ym.append(Yref)

    Xm = np.asarray(Xm)
    Ym = np.asarray(Ym)
  #  print(Xm.shape)
    Xmoy_ref,Xstd_ref=statP.average(Xref,axis=())
    
    #idea : compute the flow average over several realizations (ensemble avg)  
    Ct  = []
    for tc in range(tmin,tmax):
        Sp = []
        Xt_moy = []
        XDt_moy = []
        for M in Mlist:
            Xref,Yref=chose_axe(M,t,axes)
            Xt_m,Xstd_ref=statP.average(Xref)            
            Xt_moy.append(Xt_m)
            
            X,Y=chose_axe(M,tc,axes)
            XDt_m,Xstd=statP.average(X)
            XDt_moy.append(XDt_m)
        
        Xt_ref = statP.average(Xt_moy)
        XDt_ref = statP.average(XDt_moy)
            
        for M in Mlist:
            Xref,Yref=chose_axe(M,t,axes)
          #  Xmoy_ref,Xstd_ref=statP.average(Xref)
            X,Y=chose_axe(M,tc,axes)
          #  Xmoy,Xstd=statP.average(X)
            Sp.append((X-XDt_ref)**p*(Xref-XDt_ref)**p)
     #   print(np.shape(Sp))
        Ct.append(statP.average(Sp)[0])
     #   print(Ct)
    C0 = Ct[t-tmin]
    C_norm = np.asarray(Ct)/C0

    indices = np.argsort(t_c)
    t_c = t_c[indices]-M.t[t]
    C_norm = C_norm[indices]
    
    #plot distribution of Sp values ?
    figs = {}
    if display or save:
        field = axes[0]
        title = str(axes[0]) + ', '+str(axes[1])
        graphes.graph(t_c,C_norm,fignum=fignum,label=label)
        figs.update(graphes.legende('$t (s)$','$C$',title))    
        
        if save:
            name = 'fx_'+str(int(np.floor(M.fx*1000)))+'m_t_'+str(int(np.floor(M.t[t]*1000)))+'m_'
            graphes.save_figs(figs,prefix='./Stat_avg/Time_correlation/Serie/'+field+'/'+name,dpi=300,display=True,frmt='png')   
    return t_c,C_norm
    
def t_c(Mlist,t,field='Ux',N=100,p=1,display=False):
    for M in Mlist:
        tc,C_norm = corr_v_t([M],t,axes=[field,field],N=N,p=1,display=True,label='k^')
    
    tc,C_norm = corr_v_t(Mlist,t,axes=[field,field],N=N,p=1,display=True,label='ro')

   # print(C_norm[N+1]/2.)
    ind1 = np.argmin(np.abs(C_norm[N+1]/2.-C_norm[1:N]))
    ind2 = np.argmin(np.abs(C_norm[N+1]/2.-C_norm[N:]))
    tau = tc[ind2+N]-tc[ind1]

#    print(ind1)
#    print(ind2)
#    print(tau)
    return tau 
    
def compute_tc(Mlist,Corr_fun,fields,log=True):
    M = Mlist[0]

    labels = ['ro','k^']
    s = ''
    figs = {}
    
    for i,field in enumerate(fields):
        C = Corr_fun[field]
    
        t=[]
        tau=[]
        for tup in C:
            tc = tup[1]
            C_norm = tup[2]
            N = len(C_norm)/2
            ind1 = np.argmin(np.abs(C_norm[N+1]/2.-C_norm[1:N]))
            ind2 = np.argmin(np.abs(C_norm[N+1]/2.-C_norm[N:]))
            tau.append(tc[ind2+N]-tc[ind1])
            t.append(tc[N])
    
        if log:
            graphes.graphloglog(t,tau,fignum=1,label=labels[i])
                # graphes.set_axis(min(t),max(t),min(t),0.1)
        else:
            graphes.graph(t,tau,fignum=1,label=labels[i])
            graphes.set_axis(0,0.4,0,0.1)
        s = s + field + ', '
        figs.update(graphes.legende('t (s)','t_c (s)',s[:-2]))
                
    
        name = 'fx_'+str(int(np.floor(M.fx*1000)))+'m'
        graphes.save_figs(figs,prefix='./Stat_avg/Time_correlation/Overview/'+field+'/'+name,dpi=300,display=True,frmt='png')   
    
                
    
def corr_v_d(Mlist,t,axes=['Ux','Ux'],N=100,Dt=1,p=1,display=False,save=False,label='^',fignum=1):
    """
    compute the correlation function in space at a given time. Average over space (ie along d-1 dimensions of Mlist[axes[0]]),and eventually over time
    INPUT
    ------
    M : Mdata object
        must have attributes : t, fields ('Ux','Uy', ...)
        must have method shape()
    t : int
        time index
    axes : 2 elements string list
        attributes of M to be used for the correlation function
    N : int
        Number of frames before and after time t.
    p : int. default value is 1
        power of the fields C_p(Dt) = U1**p * U2**p / <U1**2p >
    display : bool. default value is false
    """
    #    Ct=[]    
    #how to average ?? -> need a mean on several realization ? or just substract the average value (in space) ?
    M = Mlist[0]
    tlist=M.t
   # print("Compute mean values")
    tmin = max(0,t-N)
    tmax = min(len(tlist)-1,t+N)
    #print(tmin)
    #print(tmax)
    t_c = np.asarray(M.t)[np.arange(tmin,tmax)]#-M.t[t]
   # print(np.mean(np.diff(t_c)))
    #max number of time step 
  #  print("Compute mean values 2")  
    #Compute the average flow
    Xm = []
    Ym = []
    for M in Mlist:
        Xref,Yref=chose_axe(M,t,axes)
        Xm.append(Xref)
        Ym.append(Yref)

    Xm = np.asarray(Xm)
    Ym = np.asarray(Ym)
  #  print(Xm.shape)
    Xmoy_ref,Xstd_ref=statP.average(Xref,axis=())
    
    #idea : compute the flow average over several realizations (ensemble avg)  
    Ct  = []
    for tc in range(tmin,tmax):
        Sp = []
        Xt_moy = []
        XDt_moy = []
        for M in Mlist:
            Xref,Yref=chose_axe(M,t,axes)
            Xt_m,Xstd_ref=statP.average(Xref)            
            Xt_moy.append(Xt_m)
            
            X,Y=chose_axe(M,tc,axes)
            XDt_m,Xstd=statP.average(X)
            XDt_moy.append(XDt_m)
        
        Xt_ref = statP.average(Xt_moy)
        XDt_ref = statP.average(XDt_moy)
            
        for M in Mlist:
            Xref,Yref=chose_axe(M,t,axes)
          #  Xmoy_ref,Xstd_ref=statP.average(Xref)
            X,Y=chose_axe(M,tc,axes)
          #  Xmoy,Xstd=statP.average(X)
            Sp.append((X-XDt_ref)**p*(Xref-XDt_ref)**p)
     #   print(np.shape(Sp))
        Ct.append(statP.average(Sp)[0])
     #   print(Ct)
    C0 = Ct[t-tmin]
    C_norm = np.asarray(Ct)/C0

    indices = np.argsort(t_c)
    t_c = t_c[indices]-M.t[t]
    C_norm = C_norm[indices]
    
    #plot distribution of Sp values ?
    figs = {}
    if display or save:
        field = axes[0]
        title = str(axes[0]) + ', '+str(axes[1])
        graphes.graph(t_c,C_norm,fignum=fignum,label=label)
        figs.update(graphes.legende('$t (s)$','$C$',title))    
        
        if save:
            name = 'fx_'+str(int(np.floor(M.fx*1000)))+'m_t_'+str(int(np.floor(M.t[t]*1000)))+'m_'
            graphes.save_figs(figs,prefix='./Stat_avg/Time_correlation/Serie/'+field+'/'+name,dpi=300,display=True,frmt='png')   
    return t_c,C_norm   

def structure_function(M,t,indices,axes=['E','E'],p=2):
    """
    Compute structure functions from a Mdata set
    INPUT
    -----
    M : Mdata object
    t : int 
        time index
    indices : list of tuple 
        list of pair of indices, used as the sample of Mdata to compute two points quantities
    axe : string. 
        Possible values are : 'E', 'xx', 'yy', 'xy'
        to be implemented :  'l' : longitudinal, 't' : transverse
    OUTPUT
    -----
    C : list of float
        list of evaluated structure functions on the positions specified by the list indices
    """
    X,Y=chose_axe(M,t,axes)    
    C=[]
    for i,j in indices.keys():  
        for k,l in zip(indices[(i,j)][0],indices[(i,j)][1]):
            C.append(X[i,j]**p-Y[k,l]**p)# for k,l in indices[(i,j)]])
 #   print(len(C))
    Cmoy,Cstd=statP.average(C)
    return C

def chose_axe(M,t,axes,Dt=1):
    """
    Chose N axis of a Mdata set
    INPUT
    -----
    M : Madata object
    t : int
        time index
    axes : string list 
        Possible values are : 'E', 'Ux', 'Uy', 'strain', 'omega'
    OUTPUT
    ----- 
    """
    data = tuple([access.get(M,ax,t,Dt=Dt) for ax in axes])
    return data
    
def stat_corr_t(M,t,Dt=20,axes=['Ux','Ux'],p=1,display=False,label='k^',fignum=0):

    t0 = M.t[t]
    tlist = range(t-Dt//2,t+Dt//2)
    
    curves = []
    for t in tlist:
        curves.append(corr_v_t([M],t,N=20,axes=axes,p=p,display=False))
    X,Y,Yerr = statP.box_average(curves,50)
    X = X[~np.isnan(X)]
    Y = Y[~np.isnan(Y)]
    Yerr = Yerr[~np.isnan(Yerr)]

    if display:
#        graphes.set_fig(1)
        graphes.errorbar(np.abs(X)/t0,Y,X*0,Yerr,fignum=fignum,label=label)
        graphes.legende('$t/u^{2m}$','$C_t$','$m=1/2$')

    name = 'Corr_'+axes[0]+'_'+axes[1]+'_'+str(t)
    filename = './Corr_functions/'+M.id.date+'/'+M.id.get_id() +'/'+ name + '.txt'
    
    keys = ['t',name]
    List_info = [np.ndarray.tolist(X),np.ndarray.tolist(Y)]
    
    rw_data.write_dictionnary(filename,keys,List_info,delimiter='\t')
 #   print(X)
 #   print(Y)
    return X,Y,Yerr    

def compute_Ct(M,tlist=[],axes=['Ux','Ux'],p=1,display=False,label='ko',fignum=1):
    display_part=False
    if tlist==[]:
        t0 = 20
        Dt = 50
        dimensions = M.shape()
        tlist = range(t0,dimensions[2]-t0,Dt)

    tau = np.zeros(len(tlist))
    tf = np.zeros(len(tlist))
        
    for i,t in enumerate(tlist):
        X,Y,Yerr = stat_corr_t(M,t,axes=axes,p=1,display=False)
        try:
            popt,pcurv = scipy.optimize.curve_fit(fitting.exp,np.abs(X),Y)
        except ValueError:
            print("NaN values encountered, fit skipped")
            X,Y,Yerr = stat_corr_t(M,t,axe='E',p=1,display=True)
           # input()
            pcurv=[]
            popt=[1]
        except RuntimeError:
            print("Fitting did not converge, arbitrarly chosen to previous value")
            pcurv = []
            if i==0:
                popt = [1]
            else:
                popt = [tau[i-1]]
                
        tf[i] = M.t[t]
        tau[i]=-1/popt[0]
     #   print(str(M.t[t]) + ' : ' +str(tau[i]))

        if display_part:
                texp = np.abs(X)
                graphes.set_fig(1)
                graphes.errorbar(texp,Y,texp*0,Yerr,fignum=0,label='ko')
                graphes.graph(texp,exp(texp,-1/tau[i]),fignum=1,label='r')
                graphes.legende('$t/u^{2m}$','$C_t$','$m=1/2$')
                
    if display:
        graphes.graphloglog(tf,tau,fignum=fignum,label=label)
        graphes.legende('t (s)','t_c',graphes.title(M))

    return tf,tau
    

def hist_Sp(M,t,p=2):
    """
    compute the p-th structure function <v(r)**p-v(r+d)**p>_{x,y}    
        in progress
    INPUT
    -----
        M : Mdata object to be processed
        t : int time index to process
        p : int default value : 2. Order of the structure function used for the computation
    OUTPUT
    -----
    d_para : int array
         scale of variation of the space correlation function, computed from a parabolic fit of the correlation function around 0.
         Each element correspond to a given distance d between points 
    d_lin : int array
    """
    ny,nx,nt=M.shape()
    dlist=range(1,nx)
    fignum=[4]
    label=['^-','*-','o-']
    Cxx={}
    Cyy={}
    CE={}

    indices=d_2pts(M.Ux[:,:,t],dlist[0])    
    for d in dlist:
     #   Cxx[d]=np.asarray(structure_function(M,t,indices,axe='xx',p=2))
     #   Cyy[d]=np.asarray(structure_function(M,t,indices,axe='yy',p=2))
        CE[d]=np.asarray(structure_function(M,t,indices,axe='E',p=p))

        # compute the distribution of C values 
        i=dlist.index(d)
        xbin,n=graphes.hist(np.asarray(CE[d]),Nvec=1,fignum=fignum[i],num=100,label=label[i])
        graphes.legende('$S_2 (m^2/s^2)$','pdf',str(d))
        
        nlim=100
        part_lin=np.where(np.logical_and(n>nlim,xbin>=0))
        part_para=np.where(n>nlim)

        result=np.polyfit(xbin[part_para],np.log10(n[part_para]),2,full=True)
        P_para=result[0]
        d_para=result[1]/(len(n[part_para])*np.mean(n[part_para]**2))

        result=np.polyfit(xbin[part_lin],np.log10(n[part_lin]),1,full=True)
        P_lin=result[0]
        d_lin=result[1]/(len(n[part_lin])*np.mean(n[part_lin]**2))
        
        print(d_para,d_lin)
        print("ratio = "+str(d_para/d_lin))
#        print("Curvature : "+str(P_para[0]*np.mean(xbin[part])*1000))
#        print("Straight : "+str(P_lin[0]*1000))
        CE_fit=np.polyval(P_para,xbin)
        CE_fit2=np.polyval(P_lin,xbin[part_lin])
        
        graphes.semilogy(xbin,10**CE_fit,fignum=fignum[i],label='r-')  
        graphes.semilogy(xbin[part_lin],10**CE_fit2,fignum=fignum[i],label='k--')  
        
#        print(d,len(CE[d]),CE[d])
    return d_para,d_lin5
    
def save_matrix_test():
#    time_correlation(Mlist2)
    nt = 10
    nd = 5
    t = np.arange(0,nt)
    d = np.arange(0,nd)
    
    C = np.random.rand(nd,nt)
    List_info = [d, t, C]    
    keys = ['d','t', 'Corr']
    
    filename = '/Volumes/labshared3/Stephane/Experiments/Accelerated_grid/2015_08_03/Results/Corr_test.txt'
    rw_data.write_matrix(filename,keys,List_info)
        
    #spatial_correlation(M)    

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
        

def __main__(Mlist):
    time_correlation(Mlist)

    for M in Mlist:
        spatial_correlation(M,compute=False)

#main(Mlist2)            
            