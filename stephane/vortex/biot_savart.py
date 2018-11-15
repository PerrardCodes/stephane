

import numpy as np

#import stephane.display.vfield as vfield
import stephane.display.graphes as graphes

import time
import stephane.jhtd.strain_tensor as strain_tensor
import stephane.tools.Smath as Smath
import stephane.analysis.cdata as cdata

def velocity_from_line(paths,X,Gamma=1,d=3):
    """
    Compute the 3d components velocity field at a list of specified points in space from a distribution of line vorticity
    INPUT
    -----
    paths : List of [N, 3] arrays
        coordinates of the points along the vorticity line. Each element of the list correspond to the coordinates of one closed loop
    X : List of 3 np array
        Coordinates of the point to compute the velocity field
    OUTPUT
    -----
    """
    #circulation Gamma arbitrary set to 1 for now    
    N = np.shape(X)[0]
    V = np.zeros((N,d))
    
    for path in paths:
        
        T=tangent(path)
#        print('T = ',T)
        L = np.sum(norm(T,axis=1)) #number of points : useless ?
              #  print(L/2/np.pi)
    
        for i,x in enumerate(X):
            V[i,:]=compute_u(x,Gamma,T,path,d=d)
#    print(np.min(V[:,2]))
    return V
    
    

    
def compute_u_matrix(R,Rp,Gamma,d=3):
    """
    R : locations of the point r where we compute the velocity field
    """
    
    n = Rp.shape[0]
    m = R.shape[0]

    C = np.tile(R,(n,1,1))-np.transpose(np.tile(Rp,(m,1,1)),(1,0,2))

    T = biot.tangent(Rp,d=3,step=1,cyclic=True) #compute tangent vector
    T_mat = np.transpose(np.tile(T,(m,1,1)),(1,0,2)) # make a (n,m,3) matrix
    
    dl = np.cross(T_mat,C)
    
    Cn = norm(C,axis=2)
    norm_C = np.transpose(np.asarray([Cn for k in range(d)]),(1,2,0))
    
    U = Gamma/(4*np.pi)*np.sum(dl/norm_C**3,axis=0)
    
    print(T_mat.shape)
    #print(T_mat[0,...])
    Cnorm = norm(C)
    print(C.shape)
    
    n = path.shape[0]
    m = x.shape[0]
    
    T_mat = np.tile(T,(n,1,1))
    
    C = np.tile(U,(n,1,1))-np.transpose(np.tile(U,(n,1,1)),(1,0,2))
    
    C_norm = norm(C,axis=2)
    val = np.transpose(np.asarray([C_norm for k in range(d)]))

    dl = np.cross(T,C)            
    return Gamma/(4*np.pi)*np.sum(dl/val**3,axis=0)#*L/(2*np.pi) : jacobian 
    
def compute_u(X,Gamma,T,path,d=3):
    u = x-path
    u_norm = norm(u,axis=1)
    val = np.transpose(np.asarray([u_norm for k in range(d)]))

    dl = np.cross(T,u)            
    return Gamma/(4*np.pi)*np.sum(dl/val**3,axis=0)#*L/(2*np.pi) : jacobian 
    



def norm(u,axis=1):
    return np.sqrt(np.sum(u**2,axis=axis))
    
def normalize(U):
    """
    Normalize the matrix U along its last dimension
    """
    d = len(np.shape(U))-1
    
    norm = np.sqrt(np.sum(U**2,axis=d))
    permute = tuple([i for i in range(1,d+1)]+[0])
    
    U_norm = np.transpose(np.asarray([norm for k in range(d)]),permute)
    
    return U/U_norm
    
    
def tangent(X,d=3,step=1,cyclic=True):
    """
    Compute the tangent vector on each point. 
        Assume that the first and the last point are close each other (closed loop !)
        how should we normalize it ? -> by the curviligne abscisse 
    INPUT
    -----
    X : [N,d] array
        Line of spatial coordinates
    d : int
        spatial dimension, default is 3
    step : int
        index step used to compute the spatial derivative. default is 1
    OUTPUT
    -----
    dV : 
    """
    alpha=3/4.
    beta=-3/20.
    gamma=1/60.

    N = np.shape(X)[0]
    n=3    
    
    if cyclic:
        X_ext = np.concatenate((X[N-n:,...],X,X[0:n,...]),axis=0)
        dV = np.zeros((N,d))
    else:
        X_ext = X
        dV = np.zeros((N-n*2,d))
        
    if d==1:
        X_ext = np.reshape(X_ext,(np.shape(X_ext)[0],1))
    
    #compute the first derivative along each spatial axis
    for j in range(d):
        
        tl=[[slice(3,-3,step)]+[j] for p in range(6)]
    
        tl[0][0]=slice(0,-5,step)
        tl[1][0]=slice(1,-4,step)
        tl[2][0]=slice(2,-3,step)
        tl[3][0]=slice(3,-2,step)
        tl[4][0]=slice(4,-1,step)
        tl[5][0]=slice(5,None,step)
    
        dV1 = np.sum(np.asarray([np.diff(X_ext[tl[k]],axis=0) for k in range(2,4)]),axis=0)
        dV2 = np.sum(np.asarray([np.diff(X_ext[tl[k]],axis=0) for k in range(1,5)]),axis=0)
        dV3 = np.sum(np.asarray([np.diff(X_ext[tl[k]],axis=0) for k in range(6)]),axis=0)
        
        dV[...,j] = alpha*dV1+beta*dV2+gamma*dV3
    return dV
    
def generate_vortex(R,N):
    """
    generate a circular point-like vortex
    INPUT
    -----
    R : float
        Radius of the vortex
    N : int
        number of points
    OUTPUT
    -----
    path : [N, 3] arrays
    """
    Z = 0    
    path=[]
    
    for i in range(N): 
        theta = i*2*np.pi/N
        X = R*np.cos(theta)
        Y = R*np.sin(theta)
        path.append([X,Y,Z])
        
    return np.asarray(path)
    
def B_along_axis(X,R,Gamma):
    X_norm = np.sqrt(np.sum(np.asarray(X)**2,axis=1))
  #  X_norm = np.asarray(X)[:,2]
    return Gamma*R**2/(2*(X_norm**2+R**2)**(3./2))
    
def tangent_test():
    N=100
    
    path = generate_vortex(1,N)
    
    dV = tangent(path)*N/(2*np.pi)
    
    print(np.mean(norm(dV)))
    print(np.std(norm(dV)))
    
    indices = np.arange(0,100,10)
    for i in indices:
        print(dV[i,:],path[i,:])
        
#    graphes.graph(path[:,0],path[:,1],label='r')
    graphes.graph(np.arange(N),path[:,0])
    graphes.graph(np.arange(N),np.sum(dV*path,axis=1))
#    vfield.plot(path[:,0],path[:,1],dV)
#   graphes.set_axis(-1.1,1.1,-1.5,1.5)
    
    graphes.legende('x','y','')
    
def example():
    R = 2
    Gamma = 5
    N = 10**3
    
    path = generate_vortex(R,N) #generation of a vortex ring

    eps = 0.02
    X = [np.asarray([0.,0.,z]) for z in np.arange(-10,10+eps,eps)] #axis to look at

    start=time.time()    
    V = velocity_from_line([path],X,Gamma)
    end=time.time()
    print("Time elapsed : "+str(end-start)+"s")
    print("Number of computed values :"+str(np.shape(X)[0]))

    V_th = B_along_axis(X,R,Gamma)

    x = np.asarray(X)
    
    display_profile(x,V,label='ks')
    graphes.graph(x[:,2],V_th,label='r.-',fignum=1)
    
def display_profile(x,V,label='k',axe=2,fignum=0):
    x = np.asarray(x)
    
    z = x[:,axe]
    labels = ['Ux','Uy','Uz']
    
    for i in range(3):
        graphes.graph(z,V[:,i],fignum=-i+3+fignum,label=label)
        graphes.legende(labels[i][1]+' (au)','V ',labels[i])
    
def artificial_1d():
    #just to test derivative : works perfectly
    
    dx=0.001
    p=5
    
    x = np.arange(0,1,dx)
    y = np.power(x,p)
    
    dy_num = tangent(y,d=1,step=1,cyclic=False)/dx
    dy_th = p*np.power(x[3:-3],p-1)
    
    graphes.graph(x,dy_num)
    graphes.graph(x,dy_th)
    graphes.set_axis(0,1,-1,p+1)
    
def profile(X,sigma=1,A=1):
    return A*X*np.exp(-np.power(X,2)/2/sigma**2)
 
def profile_omega(X,sigma=1,A=1):
    return A*(2-np.power(X,2)/sigma**2)*np.exp(-np.power(X,2)/2/sigma**2)
    
    
def generate(dx):
    maxv = 10
    minv = -10
    
    min2 = -10
    max2 = 10
    sigma=2.
    
    x = np.arange(minv,maxv,dx)
    y = np.arange(minv,maxv,dx)
    y = np.arange(min2,max2,dx)
    
    X,Y = np.meshgrid(x,y)
        
    R,Theta = Smath.cart2pol(X,Y)
    
    epsilon=0.01
    U = -np.sin(Theta)*profile(R,sigma=sigma)
    V = np.cos(Theta)*profile(R,sigma=sigma)
    
    dimensions = np.shape(U)+(2,)
    
    Z = np.reshape(np.transpose([V,U],(1,2,0)),dimensions)
    
    return x,y,Z
    
def artificial_vortex(dx=0.1):
    x,y,Z=generate(dx)   
    vortex_maps(x,y,Z,dx=dx,display_map=True)
    
def translate_M(M,i):
    M=cdata.rm_nan(M,'Ux',rate=0.1)
    M=cdata.rm_nan(M,'Uy',rate=0.1)
    
    M=cdata.gaussian_smooth(M,'Ux',sigma=1)
    M=cdata.gaussian_smooth(M,'Uy',sigma=1)
    
    x = M.x[0,:]
    y = M.y[:,0]
    
    Ux = M.Ux[...,i]
    Uy = M.Uy[...,i]
    
    dimensions = np.shape(Ux)+(2,)
    Z = np.reshape(np.transpose([Uy,Ux],(1,2,0)),dimensions)
    return x,y,Z
    
def vortex_maps(x,y,Z,dx=1,display_map=False):    
    dZ1 = strain_tensor.strain_tensor(Z,d=2)/dx
    
    blist = np.arange(0.5,6.,0.5)
    b0 = 2.5
    dZ2={}
    for b in blist:
        print('Compute vorticity with r='+str(b))
        dZ2[b] = strain_tensor.strain_tensor_C(Z,d=2,b=b)/dx
    
    omega1,enstrophy2 = strain_tensor.vorticity(dZ1,d=2,norm=False)
    omega2,enstrophy2 = strain_tensor.vorticity(dZ2[b0],d=2,norm=False)
    
    eigen1=strain_tensor.Lambda(dZ1,d=2)
    eigen2=strain_tensor.Lambda(dZ2[b0],d=2)
    
    strain1 = np.sqrt(np.power(eigen1['Lambda_0'],2)+np.power(eigen1['Lambda_1'],2))
    strain2 = np.sqrt(np.power(eigen2['Lambda_0'],2)+np.power(eigen2['Lambda_1'],2))
    
    figs = {}
#    dimensions = np.shape(omega1)
    
    #graphes.hist(np.reshape(omega,np.prod(dimensions)),fignum=1)
    #figs.update(graphes.legende('vorticity (s^-1)','PDF',''))
    
#    graphes.set_axis(-2.5,2.5,0,40000)
    n=3
    Xp1,Yp1 = np.meshgrid(x[n:-n],y[n:-n])
    Xp2,Yp2 = np.meshgrid(x,y)
    
    if display_map:
        for i in range(2):
            for j in range(2):
                graphes.color_plot(Xp1,Yp1,dZ1[...,i,j],fignum=i+j*2+2)
                graphes.colorbar()
                title = 'dU_'+str(i)+'/x_'+str(j)
                figs.update(graphes.legende('X','Y',title,cplot=True))

        graphes.cla(6)
        graphes.color_plot(Xp1,Yp1,-omega1,fignum=6)
        graphes.colorbar()
        figs.update(graphes.legende('X','Y','Vorticity',cplot=True))
    
        graphes.color_plot(Xp1,Yp1,strain1,fignum=7)
        graphes.colorbar()
        figs.update(graphes.legende('X','Y','Strain',cplot=True))
    
    #dissipation
        dissipation1 = np.sum(np.power(dZ1,2),axis=(2,3))
        graphes.color_plot(Xp1,Yp1,dissipation1,fignum=8)
        graphes.colorbar()
        figs.update(graphes.legende('X','Y','Dissipation',cplot=True))
    
    i1 = Xp1.shape[0]//2
    i2 = Xp2.shape[0]//2
    
    
    Xp = Xp2[i2,:]
    omega_th = profile_omega(Xp,sigma=2.)
    graphes.graph(Xp1[i1,:],omega1[i1,:],label='b-',fignum=9)
    graphes.graph(Xp2[i2,:],omega2[i2,:],label='k-',fignum=9)
    graphes.graph(Xp,omega_th,label='r--',fignum=9)    
    graphes.legende('r','vorticity','Radial vorticity profile')
    
    graphes.graph(Xp1[i1,:],strain1[i1,:],label='b.-',fignum=10)
    graphes.graph(Xp2[i2,:],strain2[i2,:],label='k^-',fignum=10)
    graphes.legende('r--','strain','Radial strain profile')
    
    j2 = Xp2.shape[1]//2
    for b in blist:
        omega2,enstrophy2 = strain_tensor.vorticity(dZ2[b],d=2,norm=False)
        graphes.graph(Xp2[i2,:],omega2[i2,:],label='k-',fignum=11)
        graphes.legende('r','vorticity','Radial vorticity profile, various b')
        
        graphes.graph([b],[omega2[i2,j2]],label='k^',fignum=12)
        graphes.legende('Circle radius r (in box unit)','max. vorticity','Influence of contour size')
        
    print(Xp.shape)
    print(Xp2.shape)
    
    #    graphes.color_plot(X,Y,R,fignum=3,vmin=0,vmax=10)
    print(figs)
    return figs
#    print(figs)
    #    print(omega)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
            
            