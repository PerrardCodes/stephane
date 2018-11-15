"""
Computation of velocity field related quantities.
Include :
- Energy 
- vorticity
- 2d divergence
- strain tensor for 2d matrices. Extension to 3d to be implemented
"""


import numpy as np
import stephane.jhtd.strain_tensor as strain_tensor
import scipy.ndimage.filters as filters
import scipy.optimize as opt

#import stephane.vortex.biot_savart as biot
#import stephane.display.graphes as graphes
import stephane.tools.Smath as Smath
import stephane.analysis.cdata as cdata

import time



def compute(M,name,step=1,Dt=1,display=False,**kwargs):
    """
    Compute a physical quantity on the Mdata
    Can be either :
    vorticity : vorticity field
    strain : strain tensor
    E : total kinetic energy
    or return any attribute of M if it already exists
    
    INPUT
    -----
    M : Mdata object
    name : string
        name of the physical quantity to be processed : 'vorticity', 'strain', 'E', or any attribute of M (Ux and Uy in particular)
    OUTPUT
    -----
    M : Mdata object
        with the field requested added to the structure
    field : string
        name of the output field ('omega' for vorticity) !!! NOT CONSISTENT : should keep the same name
    """
    
    data,field = measure(M,name,step=step,Dt=Dt,display=display,**kwargs)
    setattr(M,field,data)       
    return M,field

def measure(M,name,step=1,Dt=1,display=False,**kwargs):

    functions = {'enstrophy':enstrophy_field,'omega':vorticity_field,'vorticity':vorticity_field,'strain':strain_field,'E':get_E,'dU':dU_field}
    data=[]

#    M=cdata.gaussian_smooth(M,'Ux',sigma=0.5)    #smooth first Ux and Uy datas
#    M=cdata.gaussian_smooth(M,'Uy',sigma=0.5)

    if name in functions.keys():#hasattr(M,name):
        if name is not 'E' and display:
            print('Compute '+name)
        
        M=cdata.rm_nan(M,'Ux')    #Remove NaN values from original data
        M=cdata.rm_nan(M,'Uy')
        
        for i,t in enumerate(M.t):
            subset,field=functions[name](M,i,step=step,**kwargs)
            data.append(subset)
        data = np.asarray(data)
        dimensions = data.shape
        
     #   if len(dimensions)==3:
        tup = tuple(range(1,len(dimensions))+[0])
        #if len(dimensions)>3:        
        #tup = tuple(range(1,3)+[0]+range(3,len(dimensions))) #move the two space axis first, then time, then tensor axis    
        data = np.transpose(data,tup)
#        data = cdata.smooth(data,Dt)
    else:
        if hasattr(M,name):
            data = getattr(M,name)
            field = name
        elif 'dU_' in name:
            print('Compute tensor at a given scale')
            
            try:
                step = int(name[3:])
            except:
                print("Not a valid dU_ field")
                
            for i,t in enumerate(M.t):
                subset,field=dU_field[name](M,i,step=step,**kwargs)
                data.append(subset)
            data = np.asarray(data)          
            dimensions = data.shape
        
        #   if len(dimensions)==3:
            tup = tuple(range(1,len(dimensions))+[0])
           #if len(dimensions)>3:        
           #tup = tuple(range(1,3)+[0]+range(3,len(dimensions))) #move the two space axis first, then time, then tensor axis    
            data = np.transpose(data,tup)
        
            data = cdata.smooth(data,Dt) #smooth over time. default value is no smoothing
        else:
            print('Attribute requested unknown')
            return None,name
##depreciated            
#    if filter:
#        start = time.time()
#        M=cdata.rm_nan(M,field)
#        end = time.time()
#        print('time elapsed : '+str(end-start))
#        M=cdata.gaussian_smooth(M,field,sigma=0.5)
  #  else:
    return data,field
  #      pass
        #print('No filtering')
    
def compute_frame(M,name,frame,step=1,filter=False,Dt=1,display=False,**kwargs):
    """
    Compute a physical quantity on the Mdata, for one frame
    Can be either :
    vorticity : vorticity field
    strain : strain tensor
    E : total kinetic energy
    or return any attribute of M if it already exists
    
    INPUT
    -----
    M : Mdata object
    name : string
        name of the physical quantity to be processed : 'vorticity', 'strain', 'E', or any attribute of M (Ux and Uy in particular)
    OUTPUT
    -----
    M : Mdata object
        with the field requested added to the structure
    field : string
        name of the output field ('omega' for vorticity) !!! NOT CONSISTENT : should keep the same name
    """    
    functions = {'enstrophy':enstrophy_field,'omega':vorticity_field,'vorticity':vorticity_field,'strain':strain_field,'E':get_E,'dU':dU_field}
    data=[]
    if name in functions.keys():#hasattr(M,name):
        if name is not 'E' and display:
            print('Compute '+name)        
        data,field=functions[name](M,frame,step=step,**kwargs)
        data = cdata.smooth(data,Dt)
    else:
        if hasattr(M,name):
            data = getattr(M,name)[...,frame]
            field = name
        elif 'dU_' in name:
            print('Compute tensor at a given scale')            
            try:
                step = int(name[3:])
            except:
                print("Not a valid dU_ field")
                
            data,field=functions[name](M,frame,step=step,**kwargs)
            field = field +"_"+str(step)
            data = cdata.smooth(data,Dt)
        else:
            print('Attribute requested unknown')
            return None,name
            
    return data,field
    
     
def get_E(M,i,step=1,**kwargs):
    return np.power(M.Ux[...,i],2)+np.power(M.Uy[...,i],2),'E'
    
def strain_maps(M,i):
    x,y,Z = biot.translate_M(M,i)
    figs = biot.vortex_maps(x,y,Z)
    graphes.save_graphes(M,figs,prefix='strain_maps/',suffix='')
    
def enstrophy_field(M,i,step=1,**kwargs):
    """
    Compute the enstrophy field as the square of vorticity
    """
    omega,field= vorticity_field(M,i,step=step)
    return np.power(omega,2),'enstrophy'

def dU_field(M,i,step=1,type=2,rescale=False):
    """
    Compute the derivatives field using the stephane.jhtd.strain_tensor module
    INPUT
    -----
    M : Mdata set object
        Contain the data, either 2d or 3d in space 
    i : int
        index in time. Could be removed to compute the strain tensor both in space and time        
    OUTPUT
    -----
    omega : np array of dimension (nx-6,ny-6)
        vorticity field
    name : str
        name of the field to add to the dataset
    """
    Z,d = make_Nvec(M,i)   # Z : d+1 dimension np array 
    
    x = M.x[0,:]
    dx=np.mean(np.diff(x))
 #   print(dx)
    if dx==0 or rescale==False:
 #       print('corrected')
        dx=1
    
    if type==1:
        #print('Compute from discrete scheme in space')
        dZ = strain_tensor.strain_tensor(Z,d=d,step=step)/dx
    if type==2:
        #print('Compute from circular integral')
        dZ = strain_tensor.strain_tensor_C(Z,d=d,b=step)/dx
    
    return dZ,'dU'
#    omega,enstrophy = strain_tensor.vorticity(dZ,d=d,norm=False)
#    
#   return omega,'omega'
        
def vorticity_field(M,i,step=1,type=2,rescale=False):
    """
    Compute the vorticity field using the stephane.jhtd.strain_tensor module
    INPUT
    -----
    M : Mdata set object
        Contain the data, either 2d or 3d in space 
    i : int
        index in time. Could be removed to compute the strain tensor both in space and time        
    OUTPUT
    -----
    omega : np array of dimension (nx-6,ny-6)
        vorticity field
    name : str
        name of the field to add to the dataset
    """
    d = len(M.shape())-1
    dZ,field = dU_field(M,i,step=step,type=type,rescale=rescale)
    omega,enstrophy = strain_tensor.vorticity(dZ,d=d,norm=False)
    return omega,'omega'

def strain_field(M,i,step=1,rescale=False,type=2):
    d = len(M.shape())-1
    dZ,field = dU_field(M,i,step=step,type=type,rescale=rescale)
    eigen=strain_tensor.Lambda(dZ,d=d)
    strain = 3./d*np.sqrt(np.sum([np.power(eigen['Lambda_'+str(i)],2) for i in range(d)],axis=0))
            
    return dZ,'strain'
    
def make_Nvec(M,i):
    """
    Generate a Matrix with vector components for matrix manipulation
    """
    U = M.Ux[...,i]  #2d np array for U and V
    V = M.Uy[...,i]
  #  x = M.x[0,:]
  #  dx=np.mean(np.diff(x)) 
    d=len(np.shape(U))
    dimensions = np.shape(U)+(d,)
    #make a vector np array
    tup = tuple(range(1,d+1))+(0,)
    Z = np.reshape(np.transpose([V,U],tup),dimensions)
    
    return Z,d    #Z : vector matrix construct from U and V
    
def example(M):
    #compute vorticity for a given Mdata set
    M,name=compute(M,vorticity_field)
    M=cdata.gaussian_smooth(M,name)
    
    for i in range(200,300):
        graphes.Mplot(M,name,i)
#        input
    return M
    

    
    