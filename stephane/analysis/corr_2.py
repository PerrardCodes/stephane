


import numpy as np


import stephane.analysis.corr as corr
import stephane.analysis.statP as statP

import stephane.manager.access as access






def Time(M,field,tlist,N=50,norm_t=1.):
    """
    Compute the auto correlation function in time of field for a series of time
    Store the result in a dict format
    """
    Corr_t={}
    for i in tlist:
        t,C = corr.corr_v_t([M],i,axes=[field,field],N=N,p=1,display=False)
        Corr_t[(i,'t_'+field)]=np.asarray(t)*1./norm_t
        Corr_t[(i,'Ct_'+field)]=C    

    return Corr_t
    
def Space(M,field,tlist,N=30,Np=10**4,norm_d=1.):
    dlist = range(N)

    indices={}
    Corr_d={}

    U = access.get(M,field,0)
    
    for d in dlist:
        indices[d] = corr.d_2pts_rand(U[...,0],d,Np)

    for i in tlist:
        C=np.zeros(len(dlist))
        for d in dlist:
            C[d]=compute(M,i,indices[d],axes=[field,field])
        Corr_d[(i,'d_'+field)]=np.asarray(dlist)*1./norm_d
        Corr_d[(i,'Cd_'+field)]=C
    
    return Corr_d
    
    
def compute(M,t,indices,avg=1,axes=['U','U'],p=1,average=False):
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
        vec_t = [k-i,l-j]
        
       # Xl = project(X[i,j,...],vec_t,'l')
       # Yl = project(Y[k,l,...],vec_t,'l')
        
       # Xt = project(X[i,j,...],vec_t,'t')
       # Yt = project(Y[k,l,...],vec_t,'t')
        
        
        Sp = (X[i,j,...]-Xmoy)**p*(Y[k,l,...]-Ymoy)**p   #remove the average in space ?   -> remove by default
        C.append(Sp)# for k,l in indices[(i,j)]])  
            
        Sp_norm = ((X[i,j,...]+Y[k,l,...]-Xmoy-Ymoy)/2)**(2*p) 
        C_norm.append(Sp_norm)
            # substract the mean flow ? it shouldn't change the result so much, as there is no strong mean flow
            #-> to be checked
            # how to compute the mean flow : at which scale ? box size ? local average ?
            # Cmoy,Cstd=average(C)      
    Cf = statP.average(C)[0]/statP.average(C_norm)[0]
    return Cf

def project(U,t,typ):
    vt = t / norm(t)
    
    if typ == 'l':
        pass
    if typ == 't':
        pass
    
    
        
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