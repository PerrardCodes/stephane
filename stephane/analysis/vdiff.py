

import stephane.analysis.corr as corr
import numpy as np
import stephane.tools.Smath as Smath
import stephane.manager.access as access

import stephane.display.graphes as graphes


def Sp(Mlist,t,Dt=50,axe=['Ux','Uy'],label='k',p=1):
    
#    dlist=range(1,int(max(M.shape())/2),10)   
    dlist = [1,2,3,5,8,10,15,20,25,50]
    indices=[[] for i in range(len(dlist))]
    
    M = Mlist[0]
    dim = M.shape()
    n = np.prod(dim[:-1])
    N = 2*n#10**3
    
    #print('Compute indices lists')
    for i,d in enumerate(dlist):
        indices[i]=corr.d_2pts_rand(M.x,d,N)  
    #print('done')

    figs = {}
    for i,ind in enumerate(indices):
        print("d = "+str(dlist[i]))
        graphes.cla(i*2+1)
        graphes.cla(i*2+2)
        
        C_t = []
        C_n = []
        #print(dlist[i])
        for M in Mlist:
         #   print(M.shape())
            Ct,Cn = diff(M,t,ind,Dt=Dt,axe=axe,p=p)
                    
            graphes.distribution(Cn,label=label,normfactor=1,fignum=i*2+1,norm=False)
            figs.update(graphes.legende('C_n','PDF C_n, d'+str(dlist[i]),''))       
            
            graphes.distribution(Ct,label=label,normfactor=1,fignum=i*2+2,norm=False)
            figs.update(graphes.legende('C_t','PDF C_t, d'+str(dlist[i]),''))       
            
            C_t = C_t + Ct
            C_n = C_n + Cn
        
        n_ensemble = len(Mlist)
        
        if n_ensemble>1:
            graphes.distribution(C_n,label='r',normfactor=n_ensemble,fignum=i*2+1,norm=False)
            figs.update(graphes.legende('C_n','PDF C_n, d'+str(dlist[i]),''))       
        
            graphes.distribution(C_t,label='r',normfactor=n_ensemble,fignum=i*2+2,norm=False)
            figs.update(graphes.legende('C_t','PDF C_t, d'+str(dlist[i]),''))       
        
#        graphes.distribution(C_t,label='r',normfactor=len(Mlist),fignum=i+2,norm=False)
        
#        graphes.distribution(C_t,normfactor=1,fignum=i,norm=False)
             
    return figs
    

    
                
def diff(M,t,indices,Dt=50,avg=1,axe=['Ux','Uy'],p=1,average=False):
    #for now, only Vx and Vy gradients ... could mess up the result !!!
    #should be tangent and normal gradients in respect to d
    C=[]
    X,Y=chose_axe(M,t,axe,Dt=Dt)
    
    X = X[...,::50] #data decimation
    Y = Y[...,::50] #data decimation
        
    N = X.shape[2]
   # print(X.shape)
                            
    liste,D = get_vectors(M.x,M.y,indices)
    
    
    C_t = []
    C_n = []
        
    dx = M.x[0,1]-M.x[0,0]
    Norm = dx*np.sqrt(np.sum(np.power(D[0,:],2)))
    
    for ind,d in zip(liste,D):  
        for j in range(N):
            V1 = [X[tuple(ind[0])+(j,)],Y[tuple(ind[0])+(j,)]]
            V2 = [X[tuple(ind[1])+(j,)],Y[tuple(ind[1])+(j,)]]        
#        print(V1)
 #       print(X[ind[0]])
            Ut1,Un1 = tan_norm(V1,d) #start point
            Ut2,Un2 = tan_norm(V2,d) #stop point
        
            Sp_t = (Ut1 - Ut2)**p/Norm**p   #remove the average in space ?   -> remove by default
            Sp_n = (Un1 - Un2)**p/Norm**p   #remove the average in space ?   -> remove by default
        
            C_t.append(Sp_t)
            C_n.append(Sp_n)
         
    return C_t,C_n

def project(M,t,indices):
    dim = M.shape
    d = len(dim)
    
    tup = tuple(range(1,d-1))+(0,) #remove time from the dimensions
    
    U = np.transpose([M.Ux[...,t],M.Uy[...,t]],tup)

    U_t1 = []
    U_t2 = []
    U_n1 = []
    U_n2 = []    
    
    for ind,d in zip(liste,R):
        Ut1,Un1 = tan_norm(U[ind[0],:],d) #start point
        Ut2,Un2 = tan_norm(U[ind[1],:],d) #stop point
        
        U_t1.append(Ut1)
        U_t2.append(Ut2)
        U_n1.append(Un1)
        U_n2.append(Un2)
    
    U_t1,U_t2,U_n1,U_n2    

    U[ind,:]
    
def tan_norm(V,d):
    """
    from a maxtrix of vector U, return 
    the component tangent to d 
    the component normal to d
    
    INPUT
    -----
    U : np array with d+1 dimensions, where d is the number of space/time dimension. Last dimension length should       be 1, 2 or 3
    x,y : np arrays
        spatial coordinates
    d : 1 element dictionnary containing tuple indices
        key corresponds to the starting point, value to the ending point
    """    
    #both methods work !!
    R,Theta = Smath.cart2pol(d[0],d[1])    
    U_t1 = np.dot(V,[np.cos(Theta),np.sin(Theta)])
    U_n1 = np.dot(V,[np.sin(Theta),-np.cos(Theta)])

    t = normalize(d)    
    U_t2 = np.dot(V,t)
    U_n2 = float(np.cross(V,t))

#    print(U_t1-U_t2)
#    print(U_n1-U_n2)
    return U_t2,U_n2
    
def norm(X,p=2):
    dim = X.shape
    d = len(dim)
    Xnorm = np.sqrt(np.sum(np.power(X,p),axis=d-1))
    return Xnorm
    
def normalize(X,p=2):
        #normalize along the last axis
    Xnorm = norm(X,p=p)

    dim = X.shape
    d = len(dim)
    tup = tuple(range(1,d))+(0,)
    Xv = np.transpose(np.asarray([Xnorm for i in range(dim[d-1])]),tup)
    
    return X/Xv
    
def get_vectors(x,y,indices):
    D = []
    liste = [] 
    for start in indices.keys():
        stop = indices[start]
        liste.append([start,stop])
        D.append([x[stop] - x[start],y[stop] - y[start]])
    D = np.asarray(D)
    return liste,D
        
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