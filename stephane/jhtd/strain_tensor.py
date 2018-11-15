"""
Computation of strain tensor for arbitrary size matrices
"""


import numpy as np
import math
import stephane.tools.dict2list as dict2list

import scipy.interpolate as interpolate

import time


def strain_distribution(data,jhtd=True):
    """
    Compute eigenvectors, vorticity and geometrical quantities 
        wrap into dictionnaries the computed values from the given dataset
    INPUT 
    ----- 
    data : JHTD data format 
        data to be processed
    OUTPUT
    -----
    eigen : dict 
        contain the eigenvalues Lambda and eigenvectors lambda fields
    omega : dict 
        contain the components of the vorticity field
    cosine : dict
        contain orientation angle between lambda and omega
    """
    
    #compute the eigen values / vectors, vorticity and cosine for each time step
    omega_t = []
    eigen_t = {}
    cosine_t = {}
    
    for key in data.keys():        
        eigen,omega,cosine = geom(data[key],jhtd=jhtd)
        
        omega_t += [omega]#add_list_to_dict(omega_t,omega)
        
        eigen_t = dict2list.add_list_to_dict(eigen_t,eigen)
        cosine_t = dict2list.add_list_to_dict(cosine_t,cosine)
        
    return eigen_t,omega_t,cosine_t
    
def geom(U,d=3,jhtd=True):
    """
    Compute eigenvectors, vorticity and geometrical quantities from a d matrix U
    INPUT 
    -----
    U : d dimensions np array
        data 
    d : int
        dimension, default is 3
    OUTPUT
    -----
    eigen : dict 
        contain the eigenvalues Lambda and eigenvectors lambda fields
    omega : dict 
        contain the components of the vorticity field
    cosine : dict
        contain orientation angle between lambda and omega
    """    
    dU = strain_tensor(U,d=d,step=1,jhtd=jhtd)
    
    eigen = Lambda(dU)
    omega,enstrophy = vorticity(dU)
    W = stretching_vector(dU,omega,norm=True)
        #alignement between eigenvectors lambda_i and omega
    #    cosine = {'cosine_'+str(i): np.sum(eigen['lambda_'+str(i)]*omega,axis=3) for i in range(d)}
    cosine,eigen = alignements(eigen,omega,W,d=3)
    
    return eigen,omega,cosine
    
def alignements(eigen,omega,W,d=3):
    """
    Compute the angle defined by (lambda,omega), (lambda,W), (W, omega)
        lambda : eigenvectors of the strain tensor 
        omega : vorticity vector
        W : vortex strecthing vector
    INPUT
    -----
    eigen : dict
        contains Lambda and lambda
    omega : d+1 np array
        vorticity
    W : d+1 np array
        Vortex strecthing vector
    d : int
        dimension. default value is 3
    OUPUT
    -----
    cosine : list of dict
        cosine of the angles between (lambda,omega), (lambda,W), (W, omega)
    eigen
    """
    #alignement between eigenvectors lambda_i and vortex stretching vector
    cosine = {}
    for i in range(d):
        cosine['lambda_omega_'+str(i)] = np.sum(eigen['lambda_'+str(i)]*omega,axis=3)
        cosine['lambda_W_'+str(i)] = np.sum(W*eigen['lambda_'+str(i)],axis=3)
        cosine['W_omega']= np.sum(omega*W,axis=3)

    cosine = dict2list.dict_to_dictlist(cosine)
    eigen = dict2list.dict_to_dictlist(eigen)
    
    return cosine,eigen
    
def strain_tensor(U,d=3,step=1,jhtd=True):
    """
    Compute the strain tensor of a d dimensions matrix using a 6th order sheme in space      
    INPUT
    -----	
    U : input data, numpy array format in d+1 dimensions
        U[...,0] : x component, U[...,1] : y component, and so on
        Last dimension of U must be equal to d
        if U contains a time axis, the strain tensor calculation will also return dU_i/dt 
    d : int
        spatial dimension. default value is 3
    OUTPUT
    ------
    dU : d+2 np array
        strain tensor. The two last dimensions have shape d x d. U[...,i,j] = 
    """
	#from one time (E[:,:,:,3]), compute the values of the strain tensor. These are computed at a VERY small scale !
    #spatial dimension
    dimensions = U.shape
 #   print('Dimensions : '+str(dimensions))
    
    #6th order scheme
    n_scheme = 3*step*2
 #   print(dimensions)
    dimensions=[k-n_scheme for k in dimensions[0:d]]
    dimensions.append(d)
    dimensions.append(d)
    
    #6th order scheme
    dU = np.zeros(dimensions)

    for i in range(d):
        for j in range(d):
            dU[...,i,j]=derivative(U,i,j,step=step,d=d,jhtd=jhtd)

    return dU
    #take the symetric part
#    s_sym=(dU+np.transpose(dU,(0,1,2,4,3)))/2    
#    s_asym=(dU-np.transpose(dU,(0,1,2,4,3)))
#    print(dU[4,4,4,:,:])

def stretching_vector(dU,omega,d=3,norm=False):
    """
    Compute the stretching vector W = omega_i s_ij from the strain tensor and the vorticity
    INPUT
    -----
    dU : np array of dimensions n[0],...,n[d],d,d. (d=3 by default)
        strain tensor computed from strain_tensor function
    omega : np array of dimension n[0],...,n[d],d. (d=3 by default)
        vorticity vector
    d : int
        spatial dimension. default value is 3
    norm : bool. default False
         Normalize the result (or not)
    OUTPUT
    -----
    W : np array of dimension n[0],...,n[d],d. (d=3 by default)
        stretching vector W_i = sum((dU[i,j]omega[j]))
    """
    W=np.transpose(np.asarray([np.sum(omega*dU[:,:,:,:,j],axis=3) for j in range(d)]),(1,2,3,0))
    if norm:
        W = normalize(W)
    return W
    
def derivative(U,i,j,d=3,step=1,jhtd=False):
    """
    Compute the derivative of the j-th component of U along the i direction  
        It can be applied to a d dimensions matrix.
        A 6th order sheme in space is used to compute first spatial derivative
    INPUT
    -----	
    U : input data, numpy array format
        in any format that can be translated to numpy array by numpy.asarray
    i : int
        index of the direction of computation of derivative
    j : int
        index of the component of U to be derived
    d : int
        spatial dimension. default value is 3
    OUTPUT
    ------
    dU : numpy array. Strain tensor component s_ij of the velocity field U. 
         WARNING : a strange indexation of JHTD data forced me to introduce a switch between indices 0 and 2 for j 
    """
    U = np.asarray(U)
    #due to matrix representation of .h5 files, x and z directions are exchanged. This correction is    unfortunately dimension-dependent.        
    if jhtd and (d==3):  #in the case of JHTD, invert components.
      #  print('Permutation')
        #permutation of indices 0 and 2
        ind = [2,1,0]  
        j = ind[j]
     #   pass
        
    tl=[[slice(3,-3,step) for k in range(d)]+[j] for p in range(6)]
    
    tl[0][i]=slice(0,-5,step)
    tl[1][i]=slice(1,-4,step)
    tl[2][i]=slice(2,-3,step)
    tl[3][i]=slice(3,-2,step)
    tl[4][i]=slice(4,-1,step)
    tl[5][i]=slice(5,None,step)
    
    dU1 = np.sum(np.asarray([np.diff(U[tl[k]],axis=i) for k in range(2,4)]),axis=0)
    dU2 = np.sum(np.asarray([np.diff(U[tl[k]],axis=i) for k in range(1,5)]),axis=0)
    dU3 = np.sum(np.asarray([np.diff(U[tl[k]],axis=i) for k in range(6)]),axis=0)

    alpha=3/4.
    beta=-3/20.
    gamma=1/60.

    dU = alpha*dU1+beta*dU2+gamma*dU3
    
    #compute only the two first axis (i=0,1 and j=0,1)    
   # print(dU.shape)
    return dU
    
def strain_tensor_C(U,d=2,b=1.5,step=None):
    """
    Compute the strain tensor of a 2+1 dimension matrix using a circular integral (thanks to Dr Leonardo Gordillo)
    Value of b=1.5 /2 is the most accurate radius circle according to tests performed on experimental PIV data of vortex core size
    INPUT
    -----	
    U : input data, numpy array format in d+1 dimensions
        U[...,0] : x component, U[...,1] : y component, and so on. Last dimension of U must be equal to d
    d : int
        spatial dimension. default value is 3
    b : float
        radius of the circle used for the circular integral
    OUTPUT
    ------
    dU : d+2 np array
        strain tensor. The two last dimensions have shape d x d. U[...,i,j] = dUj / dxi
    """
	#from one time (E[:,:,:,3]), compute the values of the strain tensor. These are computed at a VERY small scale !
    #spatial dimension
    dimensions = U.shape[:-1]
    dimensions +=(d,d)
    
#    dimensions.append(d)
    
    #6th order scheme
    dU = np.zeros(dimensions)
    
    nx,ny,nd = np.shape(U)
    x = range(0,nx)
    y = range(0,ny)
    
    print(x)
    print(y)
   # print(U.shape)
    #dU = derivative_C_mat(x,y,U,d=d,b=b)

    step_by_step=True
    if step_by_step:
        fx = interpolate.RectBivariateSpline(x, y, U[...,0],kx=3,ky=3)
        fy = interpolate.RectBivariateSpline(x, y, U[...,1],kx=3,ky=3)

        for i in range(dimensions[0]):
      #  print(1.*i/dimensions[0])
            for j in range(dimensions[1]):
                val = derivative_C(U,fx,fy,i,j,d=d,b=b)
                dU[i,j,...]= val
                
        
   # print(dU.shape)
    return dU

def strain_tensor_loc(U,i,j,d=2,b=1.):
    nx,ny,nd = np.shape(U)
    x = range(0,nx)
    y = range(0,ny)
    
    fx = interpolate.RectBivariateSpline(x, y, U[...,0],kx=3,ky=3)
    fy = interpolate.RectBivariateSpline(x, y, U[...,1],kx=3,ky=3)
    
    dU = np.zeros((1,1,2,2))
    dU[0,0,...] = derivative_C(U,fx,fy,i,j,d=d,b=b)
    
    return dU
    
    
def derivative_C(U,fx,fy,x0,y0,d=2,b=1):
    """
    Alternative way of computing derivatives by using a line integral on a circle.
    Return the tensor of deformation dUij. only available in 2d for now 
    vorticity (antisymetric part of 2d tensor) and strain (symetric part of the 2d tensor)
    Could be decomposed in each component by taking the x and y component of U respectively
    INPUT
    -----
    U : numpy array of dimension 2+1
        input data, the last dimension of U is of length 2 (vector field)
    fx,fy : function with two arguments
        fx,fy should be the interpolation function of U[...,0] and U[...,1] respectively
    fx : interpolation function of U[...,0]
    OUTPUT
    -----
    dU : 2x2 numpy array
        deformation tensor dUi/dxj with (i,j) in {0,1}^2
    """
    n = int(8*b) #b should be chosen among 1/2,1,3/2, etc. so that n is divisible by 4
    ri = b
    thetai = np.arange(0,2*np.pi,2*np.pi/n)  #np.linspace(0,2*np.pi,n+1)[0:-1]
    dl = 2*np.pi*b/n  #normalisation !!!
    
    xi = x0 + ri*np.cos(thetai)
    yi = y0 + ri*np.sin(thetai)
                
    fx = np.asarray([fx(xi[i],yi[i])[0] for i in range(n)])[:,0] #evaluate the interpolated function on a circle around the point of interest
    fy = np.asarray([fy(xi[i],yi[i])[0] for i in range(n)])[:,0]

   # print(fx.shape)
    F = np.transpose(np.asarray([fx,fy]))
    Fx = np.transpose(np.asarray([fx,np.zeros(n)]))
    Fy = np.transpose(np.asarray([np.zeros(n),fy]))

    N = np.transpose(np.asarray([np.cos(thetai),np.sin(thetai)]))
    T = np.transpose(np.asarray([-np.sin(thetai),np.cos(thetai)]))

    dUx_dx = np.sum(Fx*N*dl,axis=(0,1))/(np.pi*b**2)
    dUy_dy = np.sum(Fy*N*dl,axis=(0,1))/(np.pi*b**2)
    
    dUx_dy = - np.sum(Fx*T*dl,axis=(0,1))/(np.pi*b**2)
    dUy_dx = np.sum(Fy*T*dl,axis=(0,1))/(np.pi*b**2)
    
    dU = np.zeros((2,2))
    dU[0,0] = dUx_dx 
    dU[1,1] = dUy_dy 
    dU[0,1] = dUy_dx 
    dU[1,0] = dUx_dy 
    
    return dU
    
def evaluate(X,Y,f,n=8):
    """
    Evaluate f at a arbitrary set of points of coordinates X and Y
    """
#    dim = X.shape
#    Xeval = np.ndarray.tolist(np.reshape(X,np.prod(dim)))
#    Yeval = np.ndarray.tolist(np.reshape(Y,np.prod(dim)))   
   # print(Xeval)
    feval = np.zeros((X.shape[0],Y.shape[0],n))
    for i in range(n):
        feval[...,i] = f(X[:,i],Y[:,i])
    dim=feval.shape
    feval = np.reshape(feval,tuple(dim)+(1,))
    return feval
    
def derivative_C_mat(x,y,U,d=2,b=1.):
    """
    Alternative way of computing derivatives by using a line integral on a circle.
    Return the strain tensor dUi/dxj in 2d
    INPUT
    -----
    x,y : 1d np arrays.
        x and y coordinates
    U : input data, numpy array of dimensions d+1. the last dimension of U is of length d (vector field)
    d : dimension. only implemented in 2d for now
    b : float
        circle radius used to compute the line integral. default value is 1.
    """
   # b=2.
    n = int(8*b) #b should be chosen among 1/2,1,3/2, etc. so that n is divisible by 4
    ri = b    
    dl = 2*np.pi*b/n  #normalisation !!!

    t = []
    t.append(time.time())
    fy = interpolate.RectBivariateSpline(x, y, U[...,1],kx=3,ky=3)
    fx = interpolate.RectBivariateSpline(x, y, U[...,0],kx=3,ky=3)#interpolate on a square grid. Can be evaluated only on a square grid ! (actually, compute values along a grid given two axis)
    t.append(time.time())
        
    dim = np.shape(U)
    print(dim)
    thetai = np.arange(0,2*np.pi,2*np.pi/n)  #np.linspace(0,2*np.pi,n+1)[0:-1]
    Thetai_x = np.reshape(np.tile(thetai,(dim[0],1)),(dim[0],n))
    Thetai_y = np.reshape(np.tile(thetai,(dim[1],1)),(dim[1],n))
    Thetai = np.reshape(np.tile(thetai,tuple(dim[:-1])+(1,1)),tuple(dim[:-1])+(n,1))
        
    #print(Thetai_x.shape)
    #print(Thetai_y.shape)
   # print(x.shape)
#    print(y.shape)
    
#    print(Thetai.shape)
    d = len(dim)
    X0 = np.transpose(np.tile(x,(n,1)))
    Y0 = np.transpose(np.tile(y,(n,1)))

#    print(X0.shape)
#    print(Thetai_x.shape)
#    Y0 = np.transpose(np.tile(Y,(1,n)+tuple(np.ones(2,dtype=int))),(2,3,0,1))
    Xi = X0 + ri*np.cos(Thetai_x)
    Yi = Y0 + ri*np.sin(Thetai_y)

    print(X0)
    print(Xi)
    t.append(time.time())
    
    fx = evaluate(Xi,Yi,fx,n=n)
    fy = evaluate(Xi,Yi,fy,n=n)

    t.append(time.time())
    
    dim = fx.shape
    d = len(dim)    
              
    tup = tuple(range(1,d))+(0,)          
    dim = fx.shape
    #F = np.transpose(np.asarray([fx,fy]))
    print(dim)
    Fx = np.concatenate((fx,np.zeros(dim)),axis=d-1)
    Fy = np.concatenate((np.zeros(dim),fy),axis=d-1)
    
    N = np.concatenate((np.cos(Thetai),np.sin(Thetai)),axis=d-1)
    T = np.concatenate((-np.sin(Thetai),np.cos(Thetai)),axis=d-1)
    
#    print(N.shape)
#    
#    N = np.transpose(np.asarray([np.cos(thetai),np.sin(thetai)]),tup)
#    T = np.transpose(np.asarray([-np.sin(thetai),np.cos(thetai)]),tup)
    dim = U.shape
    tup = tuple(dim[:-1])+(dim[-1],dim[-1])
    dU = np.zeros(tup)

    t.append(time.time())
    
    print(np.pi*b**2)
        
    #sum along the first two axis ?? -> probably wrong :)
    dU[...,0,0] = np.sum(Fx*N*dl,axis=(2,3))/(np.pi*b**2)
    dU[...,1,1] = np.sum(Fy*N*dl,axis=(2,3))/(np.pi*b**2)
    dU[...,1,0] = - np.sum(Fx*T*dl,axis=(2,3))/(np.pi*b**2)
    dU[...,0,1] = np.sum(Fy*T*dl,axis=(2,3))/(np.pi*b**2)
    
    
    t.append(time.time())
    
    #print(np.diff(t))
    #print(np.sum(np.diff(t)))
#    dU[0,0] = dUx_dx 
#    dU[1,1] = dUy_dy 
#    dU[0,1] = dUy_dx 
#    dU[1,0] = dUx_dy 
    
    return dU

   # print(dUy_dx-dUx_dy)    
#    print(dUx_dx+dUy_dy)
    
    #line integral :
#    t = - fx*np.sin(thetai) + fy*np.cos(thetai) 
#    n =   fx*np.cos(thetai) + fy*np.sin(thetai)    
#    omega = np.sum(t)/(np.pi*b**2)
#    div = np.sum(n)/(np.pi*b**2)
#   return dU
        

def derivative_C_GPU(U,f,x0,y0,axis,d=2,b=1.):
    """
    Alternative way of computing derivatives by using a line integral on a circle.
    step toward GPU implementation
    """
    dU = np.zeros((2,2))
    #variables to declare : n,thetai,dl,xi,yi,f,dU1,dU2
    
    n = int(8*b) #b should be chosen among 1/2,1,3/2, etc. so that n is divisible by 4
    thetai = np.arange(0,2*np.pi,2*np.pi/n)  #np.linspace(0,2*np.pi,n+1)[0:-1]
    dl = 2*np.pi*b/n  #normalisation !!!
    xi = x0 + b*np.cos(thetai)
    yi = y0 + b*np.sin(thetai)
                
    f = np.asarray([f(xi[i],yi[i])[0] for i in range(n)])[:,0] #evaluate the interpolated function on a circle around the point of interest
    if axis==0:
        dU1 = np.sum(f*np.cos(thetai)*dl,axis=(0,1))/(np.pi*b**2)
        dU2 = np.sum(-f*np.sin(thetai)*dl,axis=(0,1))/(np.pi*b**2)
        dU[0,0] = dU1
        dU[1,0] = dU2 
    if axis==2:
        dU1 = np.sum(f*np.sin(thetai)*dl,axis=(0,1))/(np.pi*b**2)
        dU2 = np.sum(f*np.cos(thetai)*dl,axis=(0,1))/(np.pi*b**2)
        dU[1,1] = dU1 
        dU[0,1] = dU2 
    
    return dU
    

def normalize(U,d=3,p=2):
    """
    Normalize a matrix U using the p-norm 
    INPUT
    -----	
    U : numpy array of dimension d+1
    d : int
        spatial dimension. default value is 3
    p : float, value for the p-norm computation. default is euclidian norm (p=2)
    OUTPUT
    -----
    U : numpy array of dimension d+1
        Normalized, so that sum(U[i,j,k,:]**2)=1
    """
    dimensions = U.shape
    d = len(dimensions)-1

    module = np.power(np.sum(np.power(U,p),axis=d),1./p)
    U = U / np.transpose(np.asarray([module for k in range(d)]),tuple(range(1,d+1)+[0]))
    
    return U
    
def Lambda(a,d=3):
    """
    Compute the eigenvalues and the eigenvector of the symetric part of a matrix a
    INPUT
    -----	
    a : numpy array of dimension d+2
        d+2 np array of dimensions (nx,ny,nz,d,d)
    d : int
        spatial dimension. default value is 3
        It might not work for other dimensions than 3 for now  
    OUTPUT
    -----
    eigen : dictionnary containing eigenvalues, eigenvectors and asymetry epsilon
        epsilon refer to the adimensionnalized value of the intermediate eigen value 
        the dictionnary contains the following fields :
        Lambda_i : (i from 0 to d) eigenvalues
        lambda_i : (i from 0 to d) eigenvectors
        epsilon
    """
    #compute the eigenvalues from the strain tensor for any point in dU2
    tup = tuple(range(d)+[d+1,d])
    s_sym=(a+np.transpose(a,tup))/2
    
    T = np.matrix.trace(s_sym,axis1=d,axis2=d+1)
    #    Tabs = np.matrix.trace(np.abs(s_sym),axis1=3,axis2=4)
    #    print(np.mean(np.mean(T/Tabs)))
        #compute eigenvalues
    N=np.prod(s_sym.shape[0:d])
    s_lin = np.reshape(s_sym,(N,d,d))
    T_lin = np.reshape(T,(N,))
    
    eigen={}
    for k in range(d):
        eigen['Lambda_'+str(k)] = []
        eigen['lambda_'+str(k)] = []        
    eigen['epsilon']=[]

    for i in range(N):
        vals,vectors = np.linalg.eigh(s_lin[i,:,:])
        vectors = vectors[:,np.argsort(vals)] 
        vals = np.sort(vals)-T_lin[i]   
        
        if d==3:
            eigen['epsilon'].append(3*vals[1]/(vals[2]-vals[0]))
        else:
            #in 2dimension, epsilon should be zero, incompressibility imply vals[0]+vals[1]=0 (if there is not out of              plane stretching !!!)
            eigen['epsilon'].append((vals[0]+vals[1])/(vals[0]-vals[1]))
                
        for k in range(d):
            eigen['Lambda_'+str(k)].append(vals[k])
            eigen['lambda_'+str(k)].append(vectors[:,k])    
            
    for key in eigen.keys():
        if key[0:6]=='Lambda':
            eigen[key]=np.reshape(np.asarray(eigen[key]),s_sym.shape[0:d]) 
        if key[0:6]=='lambda':
            eigen[key]=np.reshape(np.asarray(eigen[key]),s_sym.shape[0:d]+(d,))
        if key[0:7]=='epsilon':
            eigen[key]=np.reshape(np.asarray(eigen[key]),s_sym.shape[0:d])
        
    return eigen
    
def project(U,eigen,d=3):
    '''
    Compute the strain along the axis U, given the strain tensor expressed in a eigenbasis of decomposition
    INPUT
    -----
    U : (N,d) np array
        N vectors defining the direction of computing strain
    eigen : dict
        Contains eigenvalues and eigenvectors of a Tensor matrix (commonly the strain tensor)
    OUTPUT
    -----
    '''
    s = np.zeros(np.shape(U))
    for i in range(d):
        dotp = np.sum(eigen['lambda_'+str(i)]*U,axis=1)
        
        s[:,i]=eigen['Lambda_'+str(i)]*np.abs(dotp) #normalization factor is needed 

        mean = np.mean(s[:,i])
        print('lambda_'+str(i)+' : '+str(mean))
        
    return np.sum(s,axis=1) 
    
    
def vorticity(tau,d=3,norm=False):
    """
    Compute the vorticity from the asymetric part of the strain tensor matrix tau
    INPUT
    -----	
    tau : numpy array of dimension 3
    d : int
        spatial dimension. default value is 3
        It might not work for other dimensions than 3 for now
    OUTPUT
    -----
    vorticity : numpy array of dimensions d + 1
        the last dimension contains d values corresponding to the components of the vorticity
    enstrophy : numpy array of dimensions d
        squared modulus of the vorticity
    """    
    #compute the vorticity from the strain tensor for any point in dU2   
    tup = tuple(range(d)+[d+1,d])
    s_asym=tau-np.transpose(tau,tup)
    
    if d==3:
        vorticity = -np.asarray([s_asym[...,1,2],s_asym[...,2,0],s_asym[...,0,1]])
        shift = tuple(range(1,d+1)+[0])
        vorticity = np.transpose(vorticity,shift)
    if d==2:
        vorticity = s_asym[...,1,0]     
        
    #compute enstrophy
    dimensions = np.shape(vorticity)
    if len(dimensions)>d:
        module = np.sqrt(np.sum(np.power(vorticity,2),axis=d))
        enstrophy = np.reshape(module,dimensions[0:d])
    else:
        enstrophy = np.power(vorticity,2)
        
    if norm:     
        vorticity = normalize(vorticity,d=d,p=2)
#        dimensions = vorticity.shape
#        vorticity = dict2list.to_1d_list(vorticity,t=True)
#        vorticity = vorticity / np.transpose(np.asarray([module for k in range(d)]),(1,0))
#        vorticity = np.reshape(vorticity,dimensions)
    
    return vorticity,enstrophy
    
def divergence_2d(tau,d=2,norm=False):
    """
    Compute the 2d divergence from the symetric part of the strain tensor matrix tau
    INPUT
    -----	
    tau : numpy array of dimension 3
    d : int
        spatial dimension. default value is 3
        It might not work for other dimensions than 3 for now
    OUTPUT
    -----
    div : numpy array of dimensions d + 1
        the last dimension contains 1 value corresponding to the 2d divergence
    """    
    #compute the vorticity from the strain tensor for any point in dU2       
    if d==2:
        divergence = tau[...,0,0]+tau[...,1,1]#s_asym[...,1,0]    
        return divergence 
    else:
        print("not implemented ! (and of questionnable interest anyway for 3d data)")    
        return None
    
    
    