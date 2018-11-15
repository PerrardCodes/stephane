



"""
For various time, compute the histogramm of velocity measured, and look how it
then, should be able to launch an extrapolation in the presence of out of bound velocities

Evaluate the local gradient, or too strong spatial variation of velocity (that are not expected).
Call a spatial smoothing
"""

import numpy as np
import stephane.analysis.vgradient as vgradient
import stephane.manager.access as access

def bounds(M,t0=0):
   # t = M.t[t0]          
    if hasattr(M,'Dt'):
        ft = M.Dt[t0] # Dt is computed from indexA and indexB times ;)
    else:
        ft = M.ft
  #  print(ft)
    
#    dx = np.diff(M.x)[0,0]
#    if dx==0:
#        dx=1           
    U_noise_low = bounds_pix(M.W)[0]*M.fx/ft  # minimum value : 0.1 pixel 
    U_noise_high = bounds_pix(M.W)[1]*M.fx/ft      # maximum value : a fourth of the box size 
    
    #print(U_noise_low,U_noise_high)
#    print(M.fx)    
#    print(ft)
#    print(np.sqrt(np.nanmean(M.E[...,t0])))
#    print("")
    
    return U_noise_low,U_noise_high
    
def bounds_pix(W):
    U_low = 0.1
    U_high = W/4
    
    return U_low,U_high
    
def shear_limit(W):
    """
    Test the shear criterion : dU/W < Sopt. Sopt depends on the box size
    """
  #  frames = get_frames(M,frames)
  #  values = access.get(M,'dU',frames[0],Dt=len(frames))
  #  N = np.prod(values.shape) 
    Sopt = {16:0.3,32:0.2,64:0.05,128:0.05}  # cf Analysis and treatment of errors due to high velocity gradients in particle image velocimetry, P Meunier and T. #Leweke, Exp. in Fluids (2003)
    dUmin = -Sopt[W]*5  #0.3 for W=16, 0.2 for  W = 32, 0.05 for W = 64
    dUmax = Sopt[W]*5

    return dUmin,dUmax
    
def shear_optimum(W):
    dU_opt = shear_limit(W)[1]/5
    return dU_opt

    #convert the data in 2d
    # compute local gradients : point by point differences 

def velocity(M,frame,scale=True,display=True,W=None):
    """
    Ux = access.get(M,'Ux',frame)    
    Uy = access.get(M,'Uy',frame)    
    U = np.asarray([Ux,Uy])
    dim = U.shape
    N = np.prod(dim[1:])
    d = len(dim)
    U = np.transpose(U,tuple(range(1,d))+(0,))    
    U = np.sqrt(np.sum(np.power(U,2),axis=d-1)) #velocity modulus
    """
    U = np.sqrt(access.get(M,'E',frame))
    N = np.prod(U.shape)
    
    if scale:
        Umin,Umax = bounds(M)
    else:
        Umin,Umax = bounds_pix(W)
       
    if N==0:
        #print(U.shape)
        N = np.prod(U.shape[:-1]) 
    
    r = len(np.where(np.logical_and(U>Umin,U<Umax))[0])*100./N
    
    if display:
        print("Percentage of good values (velocity test) : "+str(r)+" %")
    return r
    
def gradient(M,frame,W=32,scale=True,display=True):
    dU = access.get(M,'dU',frame)
    
    N = np.prod(dU.shape)
    if N==0:
        #print(dU.shape)
        N = np.prod(dU.shape[:-1])
    
    
    dUmin,dUmax = shear_limit(W)
    r = len(np.where(np.logical_and(dU>dUmin,dU<dUmax))[0])*100./N

    dU_opt = shear_optimum(W)
    dU_moy = np.nanstd(dU)
    ropt = dU_moy/dU_opt
    if display:
        print("Percentage of good values (gradient test) : "+str(r)+" %")
        print("ratio measured shear / optimal value : "+str(ropt)) #greater than 1 start to be bad
    return r,ropt
 
    
def Test_v(M,frames=None,display=True,**kwargs):
    frames = get_frames(M,frames)
    r = 0.
    for frame in frames:
        r += velocity(M,frame,display=display,**kwargs)
    R = r/len(frames)
    
    if display:
        print("Percentage of good values (velocity test) : "+str(R)+" %")
        
    return R
    
def Test_dv(M,frames=None,W=32,display=True,scale=True,type=1,**kwargs):
    frames = get_frames(M,frames)
    r = 0.
    ropt = 0.
    dU = access.get(M,'dU',frames[0],Dt=len(frames),compute=False,rescale=False,type=type)
    
    for frame in frames:
        r0,ropt0 = gradient(M,frame,display=False,W=W,scale=scale)
        r += r0
        ropt +=ropt0
        
    R = r/len(frames)
    Ropt = ropt/len(frames)
    
    if display:
        import stephane.display.graphes as graphes
        dUmin,dUmax = shear_limit(W)
        
        xbin,n = graphes.hist(dU,normalize=False,num=200,range=(-0.5,0.5),**kwargs)#xfactor = Dt
        maxn = max(n)*1.2
        graphes.graph([dUmin,dUmin],[0,maxn],label='r-',**kwargs)
        graphes.graph([dUmax,dUmax],[0,maxn],label='r-',**kwargs)
        graphes.legende('','','')
    
    
        print("Percentage of good values (gradient test) : "+str(R)+" %")
        print("ratio measured shear / optimal value : "+str(Ropt)) #greater than 1 start to be bad
    
    return R   
    
def accuracy(var,M,frames=None,**kwargs):
    frames = get_frames(M,frames)

    n = min(len(frames),10)
    Ux = access.get(M,'Ux',frames[0],Dt=n)
    Uy = access.get(M,'Uy',frames[0],Dt=n)
    U = np.sqrt(access.get(M,'E',frames[0],Dt=n))
    
    dim = Ux.shape
    d = len(dim)
    
    
def v_accuracy(M,frames=None,display=True,**kwargs):
    frames = get_frames(M,frames)

    n = min(len(frames),10)
    Ux = access.get(M,'Ux',frames[0],Dt=n)
    Uy = access.get(M,'Uy',frames[0],Dt=n)
    U = np.sqrt(access.get(M,'E',frames[0],Dt=n))
    
    dim = Ux.shape
    d = len(dim)
    
    std_moy_Ux,std_Ux = compare(Ux,U,n=10,b=3)
    std_moy_Uy,std_Uy = compare(Uy,U,n=10,b=3)
        
    std_moy_U = (std_moy_Ux + std_moy_Ux)/2
    
    if display==True:
        print('Relative error velocity : '+str(std_moy_U*100)+ " %")
#    print('Relative error along y : '+str(std_moy_Uy*100)+ " %")
    
    return std_moy_U,std_Ux,std_Uy       
        
def dv_accuracy(M,frames=None,display=True,**kwargs):
    frames = get_frames(M,frames)
    n = min(len(frames),10)
    dU = access.get(M,'dU',frames[0],Dt=n)
    
    dim = dU.shape
    d = len(dim)
    
    dU2 = np.sqrt(np.sum(np.power(dU,2),axis=(d-3,d-2)))
    #dU2 = np.reshape(np.tile(dU2,dim),dim+(n,))
    
    d = len(dU2.shape)
    
#    std_moy_U1 = compare(dU[...,0,0,:],n=10,b=3)  #need to normalizeby dU2   
    std_U1 = np.nanmean([np.nanstd(dU[...,0,0,slice(i,i+3)],axis=d-1)/np.nanmean(dU2[...,slice(i,i+3)],axis=d-1) for i in range(10)]) #standard deviation along x axis
    
    std_U2 = np.nanmean([np.nanstd(dU[...,0,1,slice(i,i+3)],axis=d-1)/np.nanmean(dU2[...,slice(i,i+3)],axis=d-1) for i in range(10)])
    
    std_U3 = np.nanmean([np.nanstd(dU[...,1,0,slice(i,i+3)],axis=d-1)/np.nanmean(dU2[...,slice(i,i+3)],axis=d-1) for i in range(10)])   

    std_U4 = np.nanmean([np.nanstd(dU[...,1,1,slice(i,i+3)],axis=d-1)/np.nanmean(dU2[...,slice(i,i+3)],axis=d-1) for i in range(10)])

    std_moy_dU = np.median((std_U1+std_U2+std_U3+std_U4)/4)
    
    if display==True:
        print('Relative error velocity gradient : '+str(std_moy_dU*100)+ " %")
    
    return std_moy_dU,std_U1,std_U2,std_U3,std_U4
    
def compare(U,E,start=0,n=10,b=3):
    """
    Compare the adjacent values of U along the last axis
    """
    d = len(U.shape)
    std_U = np.nanmean([np.nanstd(U[...,slice(i,i+b)],axis=d-1)/np.nanmean(E[...,slice(i,i+b)],axis=d-1) for i in range(start,start+n)]) #standard deviation along x axis
    
    std_moy_U = np.median(std_U)
    return std_moy_U,std_U
    
def Dt_accuracy(Mdict,W=32,frames=None):
    """
    Compute the error value for adjacent values of Dt
    Mdict is a dictionnary of Mdata. The key correspond to the couple of parameters (W,Dt)
    """
    key = (1,W)
    M = Mdict[key]

    frames = get_frames(M,frames)
    n = min(len(frames),10)
        
    Dtlist = range(1,20)
    dim = Mdict[key].shape()
    
    U = np.zeros(dim[:-1]+(n,len(Dtlist)))
    Ux = np.zeros(dim[:-1]+(n,len(Dtlist)))
    Uy = np.zeros(dim[:-1]+(n,len(Dtlist)))

    dU = np.zeros(dim[:-1]+(2,2,n,len(Dtlist)))
    dU2 = np.zeros(dim[:-1]+(2,2,n,len(Dtlist)))
    
    for i,Dt in enumerate(Dtlist):
        key = (Dt,W)
        Ux[...,i] = access.get(Mdict[key],'Ux',frames[0],Dt=n)
        Uy[...,i] = access.get(Mdict[key],'Uy',frames[0],Dt=n)
        U[...,i] = np.sqrt(access.get(Mdict[key],'E',frames[0],Dt=n))
        
        dU[...,i] = access.get(Mdict[key],'dU',frames[0],Dt=n)
        
        dU_norm = np.sqrt(np.sum(np.power(dU[...,i],2),axis=(2,3)))
        dU2[...,i] = np.transpose(np.tile(dU_norm,(2,2,1,1,1)),(2,3)+(0,1,4))        
        
    for i,Dt in enumerate(Dtlist[:-1]):
        std_moy_Ux,std_Ux = compare(Ux,U,start=i,n=2,b=2)
        std_moy_Uy,std_Uy = compare(Uy,U,start=i,n=2,b=2)
        std_moy_U = (std_moy_Ux + std_moy_Ux)/2

        std_moy_dU,std_dU = compare(dU,dU2,start=i,n=2,b=2)

        print(std_moy_U,std_moy_dU)
    
    return std_moy_U
    
def W_accuracy(Mdict,Dt=1,frames=None):
    """
    Compute the error value for adjacent values of Dt
    Mdict is a dictionnary of Mdata. The key correspond to the couple of parameters (W,Dt)
    """
#    W=16
    key = (Dt,W)
    M = Mdict[key]

    frames = get_frames(M,frames)
    n = min(len(frames),10)
        
    Dtlist = range(1,20)
    dim = Mdict[key].shape()
    
    U = np.zeros(dim[:-1]+(n,len(Dtlist)))
    Ux = np.zeros(dim[:-1]+(n,len(Dtlist)))
    Uy = np.zeros(dim[:-1]+(n,len(Dtlist)))
    
    
    for i,Dt in enumerate(Dtlist):
        key = (Dt,W)
        Ux[...,i] = np.sqrt(access.get(Mdict[key],'Ux',frames[0],Dt=n))
        Uy[...,i] = np.sqrt(access.get(Mdict[key],'Uy',frames[0],Dt=n))
        U[...,i] = np.sqrt(access.get(Mdict[key],'E',frames[0],Dt=n))
    
    for i,Dt in enumerate(Dtlist[:-1]):
        std_moy_Ux,std_Ux = compare(Ux,U,start=i,n=2,b=2)
        std_moy_Uy,std_Uy = compare(Uy,U,start=i,n=2,b=2)
        
        std_moy_U = (std_moy_Ux + std_moy_Ux)/2

        print(std_moy_U)
    
    return std_moy_U
    
def get_frames(M,frames):
    if frames is None:
        return range(M.shape()[-1])
    else:
        return frames
    
    
    
    
    
    
    