

import numpy as np

#import stephane.display.vfield as vfield
import stephane.display.graphes as graphes
import stephane.display.panel as panel



def move(p,i):
    u = velocity(p,i)
    p['r'][i+1,...] = p['r'][i,...] + p['dt']*(p['noise'][i,...]+p['Gamma']*u)
    
def move_field(p,i):
    u = velocity(p,i)
    p['noise'][i,...] = get_noise(p['field'][i,...],p['r'][i,...])
    
    p['r'][i+1,...] = p['r'][i,...] + p['dt']*(p['noise'][i,...]+p['Gamma']*u)
    
def get_noise(field,r):
    
    for j in range(2):
        i  = np.argmin(r[0,j]-field[...,0])
    
def velocity(p,i):
    r = p['r'][i,:,0]-p['r'][i,:,1]
    
    u = np.cross(r/np.sum(r**2),[0,0,1])[:-1]    
#    u = np.cross(r,[0,0,1])[:-1]
    
    v = np.transpose(np.asarray([u,u]))
    
   # if i%100==0:
#        print(r)
#        print(v[:,0])
#        print(v[:,1])
    
   # print(v)
    return v
    
    
def initialize(N,epsilon=10**-2,Gamma=1,d0=1.):
    p = {}
    p['r'] = np.zeros((N,2,2))
    p['u'] = np.zeros((N,2,2))
    
    p['r'][0,:,0] = [0,d0/2.]
    p['r'][0,:,1] = [0,-d0/2.]
    
    p['dt'] = 0.2*10**-1

    p['epsilon'] = epsilon
    p['Gamma'] = Gamma
    
    p['noise'] = np.random.normal(0,p['epsilon'],(N,2,2))
    
    return p
 
def initialize_field(N,epsilon=10**-2,Gamma=1,d0=1.):
    p = {}
    p['r'] = np.zeros((N,2,2))
    p['u'] = np.zeros((N,2,2))
    
    p['r'][0,:,0] = [0,d0/2.]
    p['r'][0,:,1] = [0,-d0/2.]
    
    p['dt'] = 0.2*10**-1

    p['epsilon'] = epsilon
    p['Gamma'] = Gamma
    
    p['noise'] = np.random.normal(0,p['epsilon'],(N,2,2))
    
    return p   
    
def field():
    
    n = 100
    nt = 100
    T = 10
    B = 10
    N = 50
    
    delta = 2
    tau = 100
    
    x = np.linspace(B,-B,n)
    y = np.linspace(B,-B,n)
    t = np.linspace(0,T,nt)
    
    X,Y,T = np.meshgrid(x,y,t)
    
    R = [(R_function(B),R_function(B),R_function(B)) for i in range(N)]
    E  = np.sum([mask((X-X0)/delta,(Y-Y0)/delta,(T-T0)/tau) for X0,Y0,T0 in R],axis=0)
    
    return X,Y,T,E
    
    
def serie(Gamma,epsilon,d0=1):
    savedir = './Vortex/Advection/Gaussian_noise_linear/Stat/'
    
    n = 10**3
    N=10**4
        
    plist = []
    fignum=1
   # graphes.plt.clf()
    fig,axes = panel.make([121,122],fignum=fignum,axis='on')
    fig.set_size_inches(10,5)
    
   # print(figs)
#    print(axes)
    
    for i in range(n):
        if i%100==0:
            print(i)
        
        p = simulate(N,epsilon=epsilon,Gamma=Gamma,d0=d0)
        
        p['d'] = np.sqrt(np.sum((p['r'][...,0]-p['r'][...,1])**2,axis=1))
        p['n']=n
        p['N']=N
#        graphes.graph(p['r'][:,0,0],p['r'][:,1,0])
#        graphes.graph(p['r'][:,0,1],p['r'][:,1,1])
        
        
       # panel.sca(axes[121])
       # graphes.graph(p['r'][:,0,:],p['r'][:,1,:])
#        graphes.graph(p['r'][:,0,1],p['r'][:,1,1])

    #    panel.sca(axes[122])
    #    graphes.graph(t,p['d'])        
        
        plist.append(p)

    t = np.arange(N)*p['dt']
    
    figs = graphes.legende('X','Y','Epsilon = '+str(epsilon)+', Gamma = '+str(Gamma))
    keys = ['epsilon','Gamma','n','N']
    graphes.save_figs(figs,savedir=savedir,prefix=make_prefix(keys,p))
    
    
    R = np.asarray([p['r'] for p in plist])
    print(R.shape)
    graphes.graphloglog(t,np.std(R,axis=0)[...,1,0]**2/(epsilon**2*p['dt']),fignum=2)
    
    #graphes.graphloglog(range(N),np.std(R,axis=0)[...,1,1]/epsilon,fignum=2)

    #graphes.graphloglog(range(N),np.std(R,axis=0)[...,1,0]/epsilon,fignum=2)
    #graphes.graphloglog(range(N),np.std(R,axis=0)[...,1,1]/epsilon,fignum=2)
    #graphes.graphloglog(t,t,fignum=2,label='r--')
    
    Y = t + 2/3.*Gamma**2/d0**4*t**3#-np.arctan(t)
    graphes.graphloglog(t,Y,fignum=2,label='r--')
    
#    graphes.graphloglog([10**-1,10**1],[10**-2,10**0],fignum=2,label='r--')
#    graphes.graphloglog([10**0,10**2],[10**-1,10**5],fignum=2,label='r--')
#    graphes.graphloglog([10**3,10**5],[10**6,10**10],fignum=2,label='r--')
    
    
    figs = graphes.legende('t (#)','Variance','')
    graphes.save_figs(figs,savedir=savedir,prefix=make_prefix(keys,p))
    
    return plist,figs
    
def make_prefix(keys,p):
    s = ''
    for key in keys:
        s = s + key + '_' + str(p[key]) + '_'
    return s
    
def simulate(N,**kwargs):
    
    p = initialize(N,**kwargs)
    
    for i in range(N-1):
        move(p,i)
        
    return p
        
    
    
    
    
    
    
    
    
    
    
    
    
    
            
            