

import stephane.analysis.vgradient as vgradient
import numpy as np

def get(M,field,frame,Dt=1,Dt_filt=1,compute=False,**kwargs): 
    """
    Return a subset of the data contained in a Mdata object

    INPUT
    -----
    M : Mdata object
    field : str
        field to be extracted. Either Ux, Uy, E, omega, Enstrophy, strain
    frame : int 
        start frame index
    Dt : int
        number of frames to extract. default is 1
    Dt_filt : for filtering purpose
    compute : bool
        compute the values if True, try to recover the previoulsy stores values if False.
        default value is False
    **kwargs : key word arguments
        to be passed to vgradient.compute (computationnal option for selection the method used for instance)
    OUTPUT
    -----
    data : np array
        extracted data
    """       
    if field == 'U':
        #return both component in a vectorial format
        Ux = getattr(M,'Ux')[...,frame:frame+Dt]
        Uy = getattr(M,'Uy')[...,frame:frame+Dt]
        
        d = len(M.shape())
        tup = tuple(range(1,d+1)+[0])
        
        data = np.transpose(np.asarray([Ux,Uy]),tup)
        return data
        
    if (not hasattr(M,field)) or (compute):
        if Dt_filt>1:
            print('Filtering of the data : irreversible')
        if Dt==1:
            data_part,field = vgradient.compute_frame(M,field,frame)
            return data_part
        else:
            vgradient.compute(M,field,Dt=Dt_filt,**kwargs)
                
            data = getattr(M,field)
            data_part = data[...,frame:frame+Dt]
            return data_part
    else:
        data = getattr(M,field)
        data_part = data[...,frame:frame+Dt]
        return data_part                
    
def get_attr(M,name):
    if hasattr(M,name):
        return getattr(M,name)
    else:
        print('Attribute not found')
        return None
    
def get_all(M,field):
    dimensions = M.shape()
    Dt = dimensions[-1]
    return get(M,field,0,Dt=Dt)
    
def get_cut(M,field,s_limits,frame,Dt=1):    
    data = get(M,field,frame,Dt=Dt)
    
    indices = [slice(lim[0],lim[1]) for lim in s_limits]#+[slice(frame,frame+Dt)]
    data_cut = data[indices]
    
    return data_cut
    
def avg_vs_t(M,field,frames,Dt=10,p=2,fun=np.nanmedian):
    data = get_all(M,field)
    
    t_moy = np.zeros(len(frames))
    data_moy = np.zeros(len(frames))
    
    for i,frame in enumerate(frames):
        data_moy[i] = np.sqrt(fun(np.power(data[...,frame:frame+Dt],p)))  #compute the squared mean value
        t_moy[i] = np.mean(M.t[frame:frame+Dt])
        
    return t_moy,data_moy
    
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
    data = tuple([get(M,ax,t,Dt=Dt) for ax in axes])
    return data