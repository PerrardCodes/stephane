

import numpy as np
import tangle
import os.path

import stephane.display.graphes as graphes
import stephane.display.panel as panel

import stephane.vortex.biot_savart as biot


def noise(base,sigma,n,T=0,N=None,remove_mean=True,recompute=True,An=None):
    """
    Make a random tangle of N points with n modes   
    INPUT
    -----
    base :
        np array of dimension [N,3]
    sigma :
        float. variance of the noise
    n : 
        int. number of modes
    T:
        int. default 0. offset for the sum
    """
    
    theta,N = theta_axis(n,N=N)
    
   # An = np.ones((n,3,2))#constant coefficients
   # An = 2*sigma*np.random.random((n,3))-sigma  # uniform distribution between -sigma and sigma
    if recompute:
        An = np.random.normal(0,sigma,n*3) # gaussian distribution with width sigma
        An = np.reshape(An,(n,3))
        
        Bn = np.random.normal(0,sigma,n*3) # gaussian distribution with width sigma
        Bn = np.reshape(An,(n,3))
            
    #+An[i,j,1]*np.power(np.sin((i+T)*t),1)
    path = np.transpose(np.asarray([np.sum([base[:,j]+An[i,j]*np.cos((i+T)*theta)+Bn[i,j]*np.sin((i+T)*theta) for i in range(n)],axis=0) for j in range(3)]))
    
    if remove_mean:
        C = np.nanmean(path,axis=0)
        path = path - C
    
        norm = np.sqrt(np.sum(np.power(path,2))/N) # normalize
        path = path / norm
    
    t = tangle.Tangle(path)
   # radial_density(t)
    return t,An
    

def helicity(t):   
    t.parallel_transport_framing()
    twist = t.per_segment_twist()

    t.get_reliable_helicity()
#    t.find_crossings(t.path[0])
#    s,h = t.per_segment_stretch_helicity()
    return twist

def theta_axis(n,N=None):
    if N == None:
        dt = 1/(n*100.)
    else:
        dt = 2*np.pi/N
    theta = np.arange(0,2*np.pi,dt)
    N = len(theta)

    return theta,N

def example(p,sigma):
    n=20
    T = 0
    r0 = 1
    
    theta,N = theta_axis(n,N=None)
    base = np.asarray([[r0*np.cos(k),r0*np.sin(k),0] for k in theta])   
    t = noise(base,sigma,20)
    
    savename = './Random_path/Tests/'+str(n)+'_r0'+str(r0)+'_'+str(p+1)
    save(t,prefix=savename)

def distribution(sigma,n,display=False):
    n_p = 1
    theta,N = theta_axis(n,N=None)
    r0 = 1
    
    base = np.asarray([[r0*np.cos(k),r0*np.sin(k),0] for k in theta])   

    paths = []
    for p in range(n_p):
        #print(p)
        t = noise(base,sigma,n)
        paths.append(t.paths[0])
        
      #  h = helicity(t)
     #   graphes.hist(h,fignum=2)
        
        if p<3:
            savename = './Random_path/Tests/Examples/sigma_'+str(round(sigma*1000))+'m_n'+str(n)+'_'+str(p+1)
            save(t,prefix=savename)

    t_tot = tangle.Tangle(paths)
    if display:
        figs = radial_density(t_tot)
        figs.update(graphes.legende('R','PDF(R)',''))
 #   graphes.save_figs(figs,prefix='Random_path/Tests/R_Distributions/',suffix='_sigma_'+str(round(sigma*1000))+'m',dpi=300,display=True,frmt='png')
    return t

def radial_density(t,fignum=1,label=''):
    figs = {}
    nt = len(t.paths)    
    R_tot = []
    for j in range(nt):
        R = np.sum([t.paths[j][...,i]**2 for i in range(3)],axis=0)
        R_tot = R_tot + np.ndarray.tolist(R)
        
    graphes.hist(R_tot,log=True,fignum=fignum,label=label)
    figs.update(graphes.legende('R','PDF(R)',''))
    
    return figs    

def save(t,prefix=''):
    fn =prefix+'.tangle'
    Dir = os.path.dirname(fn)
    if not os.path.isdir(Dir):
        os.makedirs(Dir)
    t.save(fn)
    
    
def V_distribution():
    V_tot = []
    box = np.arange(-1,1,0.5)
    boxz = np.arange(-2,2,0.01)
    for i in range(1):
        t = distribution(10,20)
        for x0 in box:
            for y0 in box:
                X = [np.asarray([x0,y0,z]) for z in boxz]
                V = biot.velocity_from_line(t.paths,X,Gamma=1,d=3)
                V_tot = V_tot + np.ndarray.tolist(V)
                biot.display_profile(X,V,label='',fignum=0)
    graphes.hist(V_tot,fignum=1)
        
    
def V_distribution2(n,sigma,fignum=1):

    N = 10**4
    N_p = 10    
   # names = {name:eval(name) for name in ['sigma','n','N_p']}# this piece of code does not work systematically
    names = {'sigma':sigma,'n':n,'Np':N_p}#name:eval(name) for name in ['sigma','n','N_p']}# this piece of code does not work systematically
    
    B = 1.    
    V_tot = []
    d=3
    
    subplot = ['131','132','133']
    labels = ['Ux','Uy','Uz']
    
    fig,axes = panel.make(subplot,fignum=fignum+1,axis='on')
    fig.set_size_inches(20,5)
    figs = {}
    title = 'sigma='+str(sigma) +', n='+str(n)+', Np='+str(N_p)
    V_tot = np.zeros((N,N_p,d))
    
    for p in range(N_p):
        Theta = 2*np.pi*np.random.random(N)  # uniform distribution between 0 and 2pi
        Phi = np.pi*np.random.random(N)  # uniform distribution between 0 and pi
        R = B*np.random.random(N)
    
        X = R*np.cos(Theta)*np.sin(Phi)
        Y = R*np.sin(Theta)*np.sin(Phi)
        Z = R*np.cos(Phi)
    
        Pos = [np.asarray([x0,y0,z0]) for x0,y0,z0 in zip(X,Y,Z)]
        
        t = distribution(sigma,n)
        V = biot.velocity_from_line(t.paths,Pos,Gamma=1,d=3)
        if p==0:
            t_tot = t
        else:
            t_tot.paths = t_tot.paths + t.paths   
        
        V_tot[:,p,:] = V
    
        figs = radial_density(t,fignum=fignum,label='k')
        figs.update(graphes.legende('R','PDF(R)',title))
        
        for i,num in enumerate(subplot):
            panel.sca(axes[num])
            graphes.distribution(V[...,i],normfactor=1,a=10.,label='k',fignum=fignum+1,norm=True)
            figs.update(graphes.legende(labels[i],'PDF '+labels[i],title))
    
    figs = radial_density(t_tot,fignum=fignum,label='r')
    figs.update(graphes.legende('R','PDF(R)',title))

    for i,num in enumerate(subplot):
        panel.sca(axes[num])
        graphes.distribution(V_tot[...,i],normfactor=N_p,a=10.,label='r',fignum=fignum+1,norm=True)
        figs.update(graphes.legende(labels[i],'PDF '+labels[i],title))
    
    graphes.save_figs(figs,prefix='Random_path/V_Distributions/Compilation/',suffix=suf(names),dpi=300,display=True,frmt='png')

def vfield(sigma,modes,fignum=1):
    
    names = {'sigma':sigma,'n':modes}#name:eval(name) for name in ['sigma','n','N_p']}# this piece of code does not work 
    
    t = distribution(sigma,modes)
     
    B = 1.
    n = 100
    dx = 2*B / n
    x = np.arange(-B,B,dx)
    y = 0.
    z = np.arange(-B,B,dx)
    P = np.asarray(np.meshgrid(x,y,z))
    
    dim = P.shape
    N = np.prod(dim[1:])
    X = np.reshape(P[0,...],(N,))
    Y = np.reshape(P[1,...],(N,))
    Z = np.reshape(P[2,...],(N,))
    
    Pos = [np.asarray([x0,y0,z0]) for x0,y0,z0 in zip(X,Y,Z)]
    
    V = biot.velocity_from_line(t.paths,Pos,Gamma=1,d=3)
    
    V_2d = np.reshape(V,(n,n,3))
    E = np.sum(np.power(V_2d,2),axis=2)
    
    subplot = [111]
    fig,axes = panel.make(subplot,fignum=fignum)
    fig.set_size_inches(10,10)
    
    graphes.color_plot(x,z,E,fignum=fignum,vmin=0,vmax=0,log=True)
    c=graphes.colorbar()
    
    figs = {}
    figs.update(graphes.legende('X','Z','E'))
    
    graphes.save_figs(figs,prefix='Random_path/V_Distributions/Fields/',suffix=suf(names),dpi=300,display=True,frmt='png')
    
  #  print(np.shape(V))
    
    

def suf(args,no_points=True):
    s = ''
    greek = dict((zip([-6,-3,-2,-1,0,2,3,6],['mu','m','c','d','','h','k','M'])))
    
    for key in args.keys():
        if type(args[key])==int:
            add = str(args[key])            
        if type(args[key])==float:
            if no_points:
                ndigit = np.ceil(np.log10(args[key]))
                ndigit = ndigit - 2 + (1 if np.sign(ndigit)>0 else 0)
                i = np.argmin(np.abs(np.asarray(greek.keys())-ndigit))
                n = greek.keys()[i]
                name = greek[n]
                add = str(int(round(args[key]*10**(-n))))+name
                        
        s = s + key + '_' + add + '_'
        
    s = s[:-1]
    
 #   print(s)
    return s
    



