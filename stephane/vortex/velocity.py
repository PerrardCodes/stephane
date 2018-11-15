

from ilpm import path
import numpy as np
import stephane.vortex.biot_savart as biot
import stephane.jhtd.geom as geom
import stephane.jhtd.strain_tensor as strain

import stephane.display.graphes as graphes

def compute_velocity(grid,T,d=3):
    #Vtot = [np.zeros(3) for i in X]
    dimensions = grid.shape
    grid_1d = np.reshape(grid,(np.prod(dimensions[:-1]),3))
    
    N = grid_1d.shape[0]
    X = [grid_1d[i,:] for i in range(N)]
    
    V_grid=np.zeros(dimensions)
    for p in T.paths:
        patharr = p.path
        V = biot.velocity_from_line([patharr],X)
      
        V_grid += np.reshape(V,tuple(dimensions[:-1]+(d,)))   
        #print(V_grid.shape)
        
    return V_grid
        
def example_vortex():
    pass
    
def tangle_vs_t(fileList):
#    fileList = glob.glob(Directory+'*.tangle')
    print('Number of fileList : '+str(len(fileList)))
    
    indice = np.arange(0,700,1)
    for i,file in enumerate(fileList):
        if i in indice:
            print(i)
            T = path.load_tangle(file)
            L = T.total_length()
            
            L,epsilon = example_tangle(T,n=7,d=3)
        
            graphes.graph([indice[i]],[L],fignum=1,label='ko')
            graphes.graph([indice[i]],[epsilon],fignum=2,label='r^')
        
def example_tangle(T,n=10,d=3):
#    fn = '/Volumes/labshared/Stephane_lab1/Vortex_sim/helix_5_a16_r0bar102_6bundle_0twist_a0p08_r15_scale_100/good_traces/00000000_000.tangle'
#    T = path.load_tangle(fn)    
    patharr = T.paths[0].path

    Max = np.asarray([np.max(patharr[:,i]) for i in range(d)])
    Min = np.asarray([np.min(patharr[:,i]) for i in range(d)])
    Mean = (Max + Min)/2
    W = Max - Min
    
    D=0.01
    U = [0,100,0]
    X={}
    step=[]
    minX=[]
    maxX=[]
    for i in range(d):
        minX.append(Mean[i] - D*W[i]+U[i])
        maxX.append(Mean[i] + D*W[i]+U[i])
        step.append((maxX[i]-minX[i])/n)
        
    for i in range(d):     
        X[i] = np.arange(minX[i],maxX[i],min(step))
    
    grid = np.meshgrid(X[0],X[1],X[2])
    grid = np.transpose(grid,(1,2,3,0))

    V = np.asarray(compute_velocity(grid,T))    
    
    N = np.prod(np.shape(V)[0:3])
  #  print('Number of points : '+str(N))
    
    eigen,omega,cosine = strain.geom(V)
    figs = geom.plots(eigen,omega,cosine,10)
    
    L = T.total_length()
    epsilon = np.mean(eigen['epsilon'])
    print('L = '+str(L))
    print('Epsilon : '+str(epsilon))    
    
    return L,epsilon
#    graphes.save_figs(figs,)

#    x = np.asarray(X)
#    for i in range(3):
#        biot.display_profile(x,Vtot,axe=i,fignum=i*3)
#    np.reshape(Vtot)

def stretching_along(T,a=1,d=3):

    for k,p in enumerate(T.paths):
        path = p.path

        ds = np.mean(biot.norm(np.diff(path,axis=0)))
        print('ds = '+str(ds))
        s = np.cumsum(biot.norm(np.diff(path,axis=0)))
        print(np.shape(s))
    
        x = np.arange(-3*ds,3*ds+ds,ds)
    #matrix = np.meshgrid(x,x,x)
        noise = ds*a*(np.random.random(path.shape)-0.5)
    
 #   print(x.shape)
 #   print(path.shape)
 #   print(noise.shape)      
        matrix = np.transpose(np.asarray(np.meshgrid(x,x,x)),(1,2,3,0))
    
        X = np.zeros((7,7,7,3))
        data={}
  #  print(X.shape)
#        imax = 952
        #Z = np.zeros((noise.shape))#coordinate of the data points in space
        Z = path+noise
        
        N=len(path)
        print(N)
        for i,t in enumerate(path[:N]):
#            if i<imax:
            for j in range(d):
                X[...,j]=matrix[...,j]+t[j]+noise[i,j]
                data[i] = np.asarray(compute_velocity(X,T))    
                #Z[i,j] = t[j]+noise[i,j]

#        print(t)    
        eigen,omega,cosine = strain.strain_distribution(data)
        Tan=biot.tangent(path[:N]) # Tan gives the direction of omega (!)
        
        S = strain.project(Tan,eigen) #compute the strain along the Tan direction from strain components expressed into the eigenbasis
        #stretching_vector :

#        dU = strain.strain_tensor(data,jhtd=False)
#        W = stretching_vector(dU,omega,d=3,norm=False)
#        Vortex stretching :        
        
        epsilon = eigen['epsilon']
  #      print(epsilon)
  #      print(np.shape(epsilon))
        graphes.graph(s,epsilon[:-1],label='k.',fignum=(k+1)*a)
        graphes.legende('curvilinear coordinate s','Strain tensor asymetry','',display=False)
        
        graphes.graph(s,S[:-1],label='ro',fignum=(k+1)*a+1)
        graphes.legende('curvilinear coordinate s','Vorticity stretching','',display=False)
        
        for i in range(d):
            graphes.graph(s,eigen['Lambda_'+str(i)][:-1],fignum=(k+2)*a)
            graphes.legende('curvilinear coordinate s','Strain eigenvalues','',display=False)
            
#    print(np.shape(X))
 #   eigen,omega,cosine = strain.geom(Vf)
   # figs = geom.plots(eigen,omega,cosine,10)
        
def asymetry():
    pass
    