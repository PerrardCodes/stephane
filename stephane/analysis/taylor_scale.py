# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 16:27:38 2015

@author: stephane
"""

#measurement of the taylor scale

import numpy as np
#import stephane.analysis.strain_tensor as strain_tensor
import stephane.display.graphes as graphes

def compute(M,i,Dt=50,display=False):
    #compute the taylor scale by averaging U dU/dx over space.
    # the derivative can be taken either in the direction of U or in the direction perpendicular to it
    #call functions from the derivative module, that compute spatial derivative accurately
    nx,ny,nt = M.shape()

    start = max(0,i-Dt//2)
    end = min(nt,i+Dt//2)

    n = end-start
    Ux = M.Ux[:,:,start:end]
    Uy = M.Uy[:,:,start:end]

    #compute the strain tensor from Ux and Uy components
    edge = 3;
    d = 2
    dU = np.zeros((nx-edge*2,ny-edge*2,d,d,n))
    
    fx = max([np.mean(np.diff(M.x)),np.mean(np.diff(M.x))])  #in mm/box
    
    for k in range(n):
        U = np.transpose(np.asarray([Ux[...,k],Uy[...,k]]),(1,2,0)) #shift the dimension to compute the strain tensor along axis 0 and 1
        dU[...,k] = fx*strain_tensor.strain_tensor(U,d=2,step=1) #strain tensor computed at the box size
        
    #naive length scale, computed from Ux dUx/dx
    index =(slice(3,-3,None),slice(3,-3,None),slice(None))
    
    E_dE = Ux[index]*dU[...,0,0,:]+Uy[index]*dU[...,1,1,:]
    E = np.power(Ux[index],2)+np.power(Uy[index],2)
    
    if display:
        graphes.hist(E_dE/np.std(E_dE),num=1000,label='ko--',fignum=1)
        graphes.hist(E/np.std(E),num=1000,label='r^-',fignum=1)   
        graphes.set_axes(-10,10,1,10**5)
        graphes.legende('E','pdf(E)','')
        
    lambda_R0 = np.mean(E)/np.std(E_dE)    
    print('')
    print(str(M.t[i]) + ' : '+str(lambda_R0))
#    input()
    
    dtheta = np.pi/100
    angles = np.arange(0,np.pi,dtheta)
    
    E_dE_l = []
    E_dE_t = []
    E_theta = []

    lambda_R_l=[]
    lambda_R_t=[]

    
    for j,theta in enumerate(angles):
        U_theta = Ux[index]*np.cos(theta)+Uy[index]*np.sin(theta)   
        
        dU_l = dU[...,0,0,:]*np.cos(theta)+dU[...,1,1,:]*np.sin(theta)
        dU_t = dU[...,1,0,:]*np.cos(theta)+dU[...,0,1,:]*np.sin(theta)  #derivative of the same component, but in the normal direction
        
        #longitudinal of U dU
        E_dE_l.append(np.std(U_theta * dU_l))
        E_dE_t.append(np.std(U_theta * dU_t))
        E_theta.append(np.mean(np.power(U_theta,2)))
        
        lambda_R_l.append(E_theta[j]/E_dE_l[j])
        lambda_R_t.append(E_theta[j]/E_dE_t[j])
        
    lambda_Rl = np.mean(np.asarray(lambda_R_l))
    lambda_Rt = np.mean(np.asarray(lambda_R_t))
    
    lambda_Rl_std = np.std(np.asarray(lambda_R_l))
    lambda_Rt_std = np.std(np.asarray(lambda_R_t))
    
    print(str(M.t[i]) + ' : '+str(lambda_Rl))        
    print(str(M.t[i]) + ' : '+str(lambda_Rt))        
        
#    graphes.graph(angles,E_dE_l,fignum=1,label='ko')
#    graphes.graph(angles,E_dE_t,fignum=1,label='r^')
        
#    lambda_R = lambda_Rl
    lambdas={}
    lambdas['l_moy']=lambda_Rl
    lambdas['t_moy']=lambda_R0
    lambdas['l_std']=lambda_Rl_std
    lambdas['t_std']=lambda_Rt_std
    
    Urms = np.sqrt(np.std(E))   #E is in mm^2/s^-2
    return lambdas,Urms
        
def taylor_scale(M,fignum=1,display=True,label='k^'):
    nx,ny,nt = M.shape()    
    t = M.t
    Dt = 20
    step = 1
    
    lambda_R = {}
    Urms = []
    t_R = []
    for i in range(Dt,nt-Dt,step):
        t_R.append(t[i])
        lambdas,U = compute(M,i,Dt=Dt)
        Urms.append(U)
        if lambda_R=={}:
            for key in lambdas.keys():
                lambda_R[key]=[lambdas[key]]
        else:
            for key in lambdas.keys():
                lambda_R[key]+=[lambdas[key]]
            
    graphes.semilogx(t_R,lambda_R['t_moy'],fignum=fignum,label=label[0]+'^')       
    graphes.semilogx(t_R,lambda_R['l_moy'],fignum=fignum,label=label[0]+'>')       
    
    graphes.graphloglog(t_R,np.asarray(Urms)*np.asarray(lambda_R['t_moy']),fignum=fignum+1,label=label)   
    graphes.graphloglog(np.asarray(Urms),np.asarray(lambda_R['t_moy']),fignum=fignum+2,label=label)       
    graphes.legende('<U''>','lambda','')