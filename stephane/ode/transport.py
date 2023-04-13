import numpy as np
import pylab as plt
import scipy.integrate as inte
import scipy.special as special
import scipy.interpolate as interp

def transport(y,t,Fr):
    vx,vz = y
    a = np.sqrt(((vx-1)**2+vz**2))
    dydt = [-Fr**2*a*(vx-1),-Fr**2*a*vz-1]
    return dydt

def transport_noise(y,t,Fr,fun):#alpha = H/l, where l is the friction length
    vx,vz = y
    a = np.sqrt(((vx-1-fun(t)[0])**2+(vz-fun(t)[1])**2))
    dydt = [-Fr**2*a*(vx-1-fun(t)[0]),-Fr**2*a*(vz-fun(t)[1])-1]
    return dydt

def transport_noise_mf(y,t,Fr,c0,c1,up):#alpha = H/l, where l is the friction length
    vx,vz = y
    a = np.sqrt(((vx-1-c0)**2+(vz-c1)**2+up**2))
    
    #if np.mod(t*100,1)==0:
    #    print(t)
    dydt = [-Fr**2*a*(vx-1-c0),-Fr**2*a*(vz-c1)-1]
    return dydt

def compute_traj(theta,v0,Fr,Fr0,*args,fun=transport,tmax=10,N=10**5,H=0.05):
    y0 = [v0*np.cos(theta),v0*np.sin(theta)]#initial velocity
    t,dt = time_axis(tmax=tmax,N=N)
    sol = inte.odeint(fun, y0, t, args=(Fr,)+tuple(args))#c0 replaced by 0    
    z0 = 1/Fr0**2
    x,z = traj(sol,[0,z0],dt)
    xm,i = find_zero(x,z)
    if t[i]==tmax:
        print('traj. does not terminate')
        
    return t[:i],x[:i]*H*Fr0**2,z[:i]*H*Fr0**2

def traj_dicho(theta,v0,Fr,Fr0,sigma,Dt,N=10**4,L=1.6,H=0.05):
    nmax = 200
    i=0
    Lf = 0
    tmin = 0
    r0 = [0,1]
    y0 = [v0*np.cos(theta),v0*np.sin(theta)]#initial velocity

#    tnoise = np.arange(0,nmax)
#    w,cfun = make_noise(tnoise,sigma,0) #∆t switching noise
    #delta = 1
    h = 1
    T,X,Z = [],[],[]
    
    while i<nmax and (r0[1]-(1-1/Fr0**2))*H*Fr0**2>0 and r0[0]<L/H/Fr0**2:
        tmax = tmin+Dt
        t,dt = time_axis(tmin=tmin,tmax=tmax,N=N)

        wp = np.random.normal(0,sigma)#*v_profile(1,1,h) #H0 = delta
        s = inte.odeint(transport_noise_mf, y0, t, args=(Fr,0,wp,0))#c0 replaced by 0    

        x,z = traj(s,r0,dt)
        r0 = [x[-1],z[-1]]
        y0 = s[-1,:]
        tmin = tmax

        T.append(t)
        X.append(x*H*Fr0**2)
        Z.append((z-(1-1/Fr0**2))*H*Fr0**2)
        
    xm,i = find_zero(X[-1],Z[-1])
    #print(len(X[-1][:i]))
    X[-1]=X[-1][:i+1]
    Z[-1]=Z[-1][:i+1]
    T[-1]=T[-1][:i+1]
        
#    T = np.asarray(T)
#    X = X*H*Fr0**2
#    Z = (Z-(1-1/Fr0**2))*H*Fr0**2
    
    return T,X,Z

def traj_beam_dicho(thetalist,vlist,Fr,Fr0,*args,display=False,N=10**5,H=0.05,L=1.6):
    Xlist,Zlist = [],[]
    X = []
    Z = []
    for j,(theta,v0) in enumerate(zip(thetalist,vlist)):
        t,x,z = traj_dicho(theta,v0,Fr,Fr0,*args,N=N)
        
        X.append(x)
        Z.append(z)
        
        if display:
            plt.plot(x,z)
            graphes.legende('$x$ (m)','$z$ (m)','')

        xmax,zL = get_measure(x,z,L=L)
        Xlist.append(xmax)
        Zlist.append(zL)

            
    Xlist = np.asarray(Xlist)
    Zlist = np.asarray(Zlist)
    return Xlist,Zlist,X,Z

def get_measure(x,z,L=1.6):
    #look for the maximum X in the last portion of trajectory
    xmax = np.max(x[-1])
    #look for the Z at X=L
    if xmax>=L:
        index = np.argmin(np.abs(x[-1]-L))
        zL = z[-1][index]
    else:
        zL = 0
    return xmax,zL

def display_beam_dicho(X,Z,xref,zref,ax=None,title=''):
    import stephane.display.graphes as graphes
    for j,(x,z) in enumerate(zip(X,Z)):
        for (xp,zp) in zip(x,z):
            ax.plot(xp,zp)
        
    for (xp,zp) in zip(xref,zref):
        ax.plot(xp,zp,'k-',linewidth=3)

    plt.vlines(1.6,0,0.2,colors='k',linestyles='dashed',alpha=1)
#    plt.hlines(0,0,3,colors='k',linestyles='solid',alpha=1)

    plt.axis([0,3,0,0.2])

    figs = graphes.legende('$x$ (m)','$z$ (m)',title,ax=ax,cplot=True)
    return figs


def v_profile(U,delta,x):
    return 0.1*U*np.exp(-x/delta)

def traj_noise_deprecated(theta,Fr,Fr0,sigma,nc,display=False,N=10**5):
    y0 = [v0*np.cos(theta),v0*np.sin(theta)]#initial velocity
    t = np.linspace(0,10,N+1)
    dt = t[1]-t[0]
    tnoise = np.linspace(0,11,int(N*1.1+1))
    
    c,cfun = make_noise_2d(tnoise,sigma,nc)

    dt = t[1]-t[0]
    sol = inte.odeint(transport_noise, y0, t,args=(cfun,Fr),hmin=dt)#c0 replaced by 0
    x,z = traj(sol,[0,1],dt)#get the trajectory from the velocity 
    xm,i = find_zero(x,z)
        
    #print(theta,t[i])

    if display:
        plt.plot(x[:i],z[:i])
    return x[:i],z[:i],xm

def find_zero(x,z):
    if np.sum(z<0)==0:
        i = np.argmin(np.abs(z))
    else:
        while np.sum(z<0)>0:
            j = np.argmin(z)
            z = z[:j]
            i = np.argmin(np.abs(z))
            z = z[:i]                
    return x[i],i

def get_height(X,Z,L=1.6):
    inds = [np.argmin(np.abs(x-L)) for x in X]
    Zlist = [z[inds[i]] for i,z in enumerate(Z)]
    return Zlist

def traj(y,r0,dt):
    R = np.cumsum(y,axis=0)*dt
    x = R[:,0]+r0[0]
    z = R[:,1]+r0[1]
    return x,z

def time_axis(tmin=0,tmax=10,N=10**5):
    t = np.linspace(tmin,tmax,N+1)
    dt = t[1]-t[0]
    return t,dt
    
def traj_beam(thetalist,vlist,Fr,Fr0,*args,fun=transport,display=False,N=10**5,H=0.05):
    Xlist = []
    X = []
    Z = []
    for j,(theta,v0) in enumerate(zip(thetalist,vlist)):
        t,x,z = compute_traj(theta,v0,Fr,Fr0,*args,fun=fun,N=N,H=H)
        
        X.append(x)
        Z.append(z)
        
        if display:
            plt.plot(x,z)
            graphes.legende('$x$ (m)','$z$ (m)','')
        Xlist.append(np.max(x))
    return np.asarray(Xlist),X,Z

def traj_beam_adv(thetalist,vlist,wlist,Fr,Fr0,u0,up,fun=transport_noise_mf,display=False,N=10**5,H=0.05,tmax=10):
    Xlist = []
    X = []
    Z = []
    for j,(theta,v0,w0) in enumerate(zip(thetalist,vlist,wlist)):
        t,x,z = compute_traj(theta,v0,Fr,Fr0,u0,w0,up,fun=fun,N=N,H=H,tmax=tmax)
        
        if len(x)==0:
            print(j,'no trajectory computed !')
            print(theta,v0,Fr,Fr0,u0,w0,up,N)
            t,x,z = compute_traj(theta,v0,Fr,Fr0,u0,w0,up,fun=fun,N=N*10,H=H,tmax=tmax*10)
            print(x)
            x = [0]
            z = [0]#Z[-1][0]
        X.append(x)
        Z.append(z)
        
        if display:
            plt.plot(x,z)
            graphes.legende('$x$ (m)','$z$ (m)','')
        Xlist.append(np.max(x))
    return np.asarray(Xlist),X,Z

def distribution(thetalist,Xlist,sigma,display=False):#starting from a distribution of angle, what is the distribution of time of flight ?
    rho = 1/np.sqrt((2*np.pi*sigma**2))*np.exp(-thetalist**2/(2*sigma**2))
    dtheta = thetalist[1]-thetalist[0]
    if display:
        plt.semilogy(thetalist,rho,'ko')

    dTh = dtheta/np.diff(Xlist)
    rhom = (rho[1:]+rho[:-1])/2
    rho_x = rhom*dTh
    xm = (Xlist[1:]+Xlist[:-1])/2
    thetam = (thetalist[1:]+thetalist[:-1])/2
    if display:
        plt.figure()
        plt.semilogy(xm,rho_x,'ro')
    return xm,rho_x

#def transport_noise_real(y,t,fun,Fr):#alpha = H/l, where l is the friction length
#    vx,vz = y
#    a = np.sqrt(((vx-Fr*(1+fun(t)[0]))**2+(vz-Fr*fun(t)[1])**2))
#    dydt = [-alpha*a*(vx-Fr*(1+fun(t)[0])),-alpha*a*(vz-Fr*fun(t)[1])-1]
#    return dydt

def make_noise(t,sigma,nc):
    n = len(t)
    a = np.random.normal(0,sigma,size=n+nc)
    if nc>0:
        b = np.cumsum(a,axis=0)
        c = b[nc:,:]-b[:-nc,:]
        c = c/np.sqrt(nc)
    else:
        c = a
    cfun = interp.interp1d(t, c, kind='cubic')
    return c,cfun

#create the noise
def make_noise_2d(t,sigma,nc):
    n = len(t)
    a = np.random.multivariate_normal([0.,0.],[[sigma**2,0],[0,sigma**2]], size=n+nc,)
    b = np.cumsum(a,axis=0)
    c = b[nc:,:]-b[:-nc,:]
    c = c/np.sqrt(nc)
    
    cf1 = interp.interp1d(t, c[:,0], kind='cubic')
    cf2 = interp.interp1d(t, c[:,1], kind='cubic')
    cfun = lambda x:(cf1(x),cf2(x))
    return c,cfun

def traj_beam_noise(thetalist,*args,**kwargs):
    Xlist = []
    for j,theta in enumerate(thetalist):
        x,z,xm = traj_noise(theta,args,kwargs)
        Xlist.append(xm)
    return np.asarray(Xlist)

def plot_trajectories(data,axs,Fr,Fr0,up,vlist,thetalist,color='w--'):
    fx = 0.132 
    H = 0.05
    Hmm = H*1000

    xc = data['x0']
    yc = data['y0']

    Xlist,X,Z = traj_beam(thetalist,vlist,Fr,Fr0,up,display=False,N=10**5,H=H,tmax=20)

    for i,k in enumerate(['near','far']):
        if k=='near':
            xf = 0
        else:
            xf = -1491/fx#distance in pix

        for x,z in zip(X,Z):
            axs[i].plot(x/fx+xc+xf,-z/fx+yc+Hmm/fx,color)
    #plt.axis([0,4000,2500,0])

    #figs = {}
    #figs.update(graphes.legende('$X$ (pixel)','$Z$ (pixel)','$U_0$ = '+str(U0)+' m/s'))
    
#graphes.save_figs(figs,savedir=savefolder,frmt='png',prefix = 'fit_U0'+str(int(key*100)))
