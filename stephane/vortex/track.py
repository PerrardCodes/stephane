
import numpy as np
import scipy.ndimage.filters as filters
import scipy.optimize as opt
import scipy.interpolate as interpolate

import stephane.manager.access as access

import stephane.vortex.core_size as core
import stephane.analysis.vgradient as vgradient
import stephane.jhtd.strain_tensor as strain

import stephane.display.graphes as graphes

import stephane.tools.browse as browse

def position(M,field='omega',indices=None,step=1,sigma=1.):
    """
    Measure the position of the min and max of vorticity over time, + their value, + the associated circulation
    Return a dictionnary containing all the values
    """

    Y = access.get_all(M,field)    
    nx,ny,nt = M.shape()
        
    if indices is None:
        indices = np.arange(step-1,nt-step+1,step)
    
    pos = {}
    variables = ['t','Xmin','Xmax','Ymin','Ymax','Vmin','Vmax','Gammamin','Gammamax']
    for var in variables:
        pos[var]=[]
    
    for i in indices:
       # print(i)
        Z = Y[:,5:,i]

        pos['t'].append(M.t[i])

        xmin,xmax,ymin,ymax = positions(M,i,sigma=sigma)    
        vmin,vmax,G = amplitude(M,i)
        
        pos['Gammamin'].append(G[0])
        pos['Gammamax'].append(G[1])

        pos['Vmin'].append(vmin)
        pos['Vmax'].append(vmax)

        pos['Xmin'].append(xmin)
        pos['Xmax'].append(xmax)
        pos['Ymin'].append(ymin)
        pos['Ymax'].append(ymax)

    return pos

def positions(M,i,field='omega',display=False,sigma=1.,**kwargs):
    """
    find the maximum position of the vortex core at a given instant, using the field specified in argument
    """

    imin,imax,jmin,jmax,Zfilt,X,Y = position_index(M,i,field=field,display=display,sigma=sigma,**kwargs)
    
    deltax_max,deltay_max = SubPix2DGauss(Zfilt/np.max(Zfilt),imax,jmax,SubPixOffset=0)
    deltax_min,deltay_min = SubPix2DGauss(Zfilt/np.min(Zfilt),imin,jmin,SubPixOffset=0)
    
 #   print(deltax_max)
 #   print(deltax_min)
 #   print(deltay_min)
    
#    print(deltax_max)
    #if deltax_max == np.nan or deltax_min==np.nan:
    #    return np.nan,np.nan,np.nan,np.nan
    #else:
    xmin = X[imin,jmin]+deltax_min
    xmax = X[imax,jmax]+deltax_max
    ymin = Y[imin,jmin]-deltay_min
    ymax = Y[imax,jmax]-deltay_max


  #  print(Zfilt.shape)
  #  print(Y.shape)
  #  print(X.shape)
    
    if display:
        graphes.color_plot(X[:,5:],Y[:,5:],Zfilt/np.min(Zfilt),**kwargs)
        graphes.graph([xmin+10],[ymin],label='bo',**kwargs)
        graphes.graph([xmax+10],[ymax],label='ro',**kwargs)
 
  #  graphes.color_plot(X[:,5:],Y[:,5:],Zfilt/np.max(Zfilt),fignum=2)
    
    #input()
#    t.append(M.t[i])
#    Xmin.append(M.x[imin,jmin])
#    Xmax.append(M.x[imax,jmax])
#    Ymin.append(M.y[imin,jmin])
#    Ymax.append(M.y[imax,jmax])
    
    return xmin,xmax,ymin,ymax
#    core.normalize(X,Y,Z)

def position_index(M,i,field='omega',display=True,sigma=1.,**kwargs):    
    dx = M.x[0,1]-M.x[0,0]
    
#    x0 = M.x[0,5]
#    y0 = M.y[0,0]
    nx,ny,nt = M.shape()
    crop=True
    if crop:
        lin = slice(0,nx)
        col = slice(0,ny)
        X = M.x[lin,col]
        Y = M.y[lin,col]
    
    Zfilt = smoothing(M,i,field=field,sigma=sigma)
    if display:
        graphes.color_plot(X,Y,Zfilt,**kwargs)

    imax,jmax = np.unravel_index(np.argmax(Zfilt),Zfilt.shape)
    imin,jmin = np.unravel_index(np.argmin(Zfilt),Zfilt.shape)
    
    return imin,imax,jmin,jmax,Zfilt,X,Y
    

def smoothing(M,i,field='omega',sigma=1.):
    Z = access.get(M,field,i)
    #Z = Z[:,5:]  #crop image
    return filters.gaussian_filter(Z, sigma=sigma)[...,0]
    
def amplitude(M,i,field='omega',display=False):
    
    xmin,xmax,ymin,ymax = positions(M,i,field='omega',display=False)

    Zfilt = smoothing(M,i,field='omega',sigma=3.)    

    x = M.x[0,:]
    y = -M.y[:,0]
    dx = M.x[0,1]-M.x[0,0]    
    Omegai = interpolate.RectBivariateSpline(y, x, Zfilt,kx=3,ky=3)
    
    U,d = vgradient.make_Nvec(M,i)   # Z : d+1 dimension np array 
    b0 = 10.
    
    imin,imax,jmin,jmax,Zfilt,X,Y = position_index(M,i,field=field,display=display)
    
    tau1 = strain.strain_tensor_loc(U,imin,jmin,d,b=b0)
    tau2 = strain.strain_tensor_loc(U,imax,jmax,d,b=b0)
    
    G = []
    for tau in [tau1,tau2]:
        omega,enstrophy = strain.vorticity(tau,d=2,norm=False)
        div = strain.divergence_2d(tau,d=2)
        G.append((omega[0,0]) * np.pi*b0**2*dx**2)  #-div[0,0]

#    divergence.append(div[0,0]/np.abs(omega[0,0]))
    val_min = Omegai(-ymin,xmin)
    val_max = Omegai(-ymax,xmax)
    
    return val_min[0,0],val_max[0,0],G
    
    
def SubPix2DGauss(c,i,j,SubPixOffset=0):
   # sub-pixel displacement using 2 dimensionnal Gaussian peak fit
    a=1
    c_p=c[i-a:i+a+1,j-a:j+a+1]
    
    if np.shape(c_p)==(3,3):#((cx.size>2) & (cy.size>2)): #(cx.size>2 & cy.size>2):
        c10=np.zeros((3,3))
        c01=np.zeros((3,3))
        c11=np.zeros((3,3))
        c20=np.zeros((3,3))
        c02=np.zeros((3,3))
        
        for i in range(-1,2):
            for j in range(-1,2):
                #following 15 lines based on
                #H. Nobach & M. Honkanen (2005)
                #Two-dimensional Gaussian regression for sub-pixel displacement
                #estimation in particle image velocimetry or particle position
                #estimation in particle tracking velocimetry
                #Experiments in Fluids (2005) 38: 511515                            
                c10[j+1,i+1]=i*np.log(c_p[i+1, j+1])
                c01[j+1,i+1]=j*np.log(c_p[i+1, j+1])
                c11[j+1,i+1]=i*j*np.log(c_p[i+1, j+1])
                c20[j+1,i+1]=(3*i**2-2)*np.log(c_p[i+1, j+1])
                c02[j+1,i+1]=(3*j**2-2)*np.log(c_p[i+1, j+1])
                
        c10=(1./6)*np.sum(c10)
        c01=(1./6)*np.sum(c01)
        c11=(1./4)*np.sum(c11)
        c20=(1./6)*np.sum(c20)
        c02=(1./6)*np.sum(c02)
        
        
        deltax=(c11*c10-2*c01*c20)/(4.*c20*c02-c11**2)
        deltay=(c11*c01-2*c10*c02)/(4.*c20*c02-c11**2)
        
        return (deltax,deltay)
    else:
        return OutOfBound()
#    peakx=x+deltax;
#    peaky=y+deltay;
    
#    SubpixelX=peakx-(interrogationarea/2)-SubPixOffset;
#    SubpixelY=peaky-(interrogationarea/2)-SubPixOffset;
#    vector=[SubpixelX, SubpixelY];

def OutOfBound():
 #   print('Index out of bound ! Keep initial value')
    return (np.nan,np.nan)
    
def profile_avg(Mlist,i,display=False,fignum=1):
    field = 'omega'
    
    R_list=np.arange(0.25,15.,0.25)
    R = np.arange(0.,15.,0.25)
    n = len(Mlist)
    
    figs = {}
    Gamma = np.zeros((n,len(R),2))
    Flux = np.zeros((n,len(R),2))
    
    for k,M in enumerate(Mlist):
        x = M.x[0,5:]
        y = -M.y[:,0]
        dx = M.x[0,1]-M.x[0,0]
        
        U,d = vgradient.make_Nvec(M,i)   # Z : d+1 dimension np array 
        
        imin,imax,jmin,jmax,Zfilt,X,Y = position_index(M,i,field=field,display=display)
        
        G = [[] for j in range(len(R_list)+1)]
        G[0].append(0)
        G[0].append(0)

        D = [[] for j in range(len(R_list)+1)]
        D[0].append(0)
        D[0].append(0)
        
        for j,b in enumerate(R_list):
            tau1 = strain.strain_tensor_loc(U,imin,jmin,d,b=b)
            tau2 = strain.strain_tensor_loc(U,imax,jmax,d,b=b)
    
            for tau in [tau1,tau2]:
                omega,enstrophy = strain.vorticity(tau,d=2,norm=False)
                div = strain.divergence_2d(tau,d=2)
                G[j+1].append((omega[0,0]) * np.pi*b**2*dx**2)  #-div[0,0]
                D[j+1].append((div[0,0]) * np.pi*b**2*dx**2)  #-div[0,0]
                
        Gamma[k,...] = np.asarray(G)
        Flux[k,...] = np.asarray(D)
        
        graphes.graph(R,-Gamma[k,:,0],label = 'kv',fignum=fignum*2)
        graphes.graph(R,Gamma[k,:,1],label = 'k^',fignum=fignum*2)
        
        graphes.graph(R,Flux[k,:,0],label = 'kv',fignum=fignum*2+1)
        graphes.graph(R,Flux[k,:,1],label = 'k^',fignum=fignum*2+1)
        
    
    Gamma_moy = np.nanmean(Gamma,axis=0)
    Flux_moy = np.nanmean(Flux,axis=0)
    
    
    graphes.graph(R,-Gamma_moy[:,0],label = 'r--',fignum=fignum*2)
    graphes.graph(R,Gamma_moy[:,1],label = 'r--',fignum=fignum*2)
    graphes.graph(R,np.nanmean(np.asarray([Gamma_moy[:,1],-Gamma_moy[:,0]]),axis=0),label = 'rs',fignum=fignum*2)
    figs.update(graphes.legende('Distance to center (mm)','Circulation (mm^2/s)',''))
    
    graphes.graph(R,Flux_moy[:,0],label = 'r--',fignum=fignum*2+1)
    graphes.graph(R,Flux_moy[:,1],label = 'r--',fignum=fignum*2+1)
    graphes.graph(R,np.nanmean(np.asarray([Flux_moy[:,1],Flux_moy[:,0]]),axis=0),label = 'rs',fignum=fignum*2+1)
    figs.update(graphes.legende('Distance to center (mm)','Divergence (mm^2/s)',''))
    
    savedir = title(Mlist[0])    
    graphes.save_figs(figs,savedir=savedir,suffix='',prefix='frame_'+str(i),frmt='pdf',dpi=300,display=True)
    
    M_profile = np.nanmean(np.asarray([Gamma_moy[:,1],-Gamma_moy[:,0]]),axis=0)
    return np.mean(M_profile[-10:]),np.std(Gamma[:,-10:,1])
    
    
def title(M):
    date = M.Id.date    
    typ = browse.get_string(M.dataDir,'piston12mm_',end='_f5Hz')
    savedir = './Vortex_Turbulence/Vortex_propagation/'+date+'/'+typ+'/'
    print(savedir)
    return savedir

def correct(u,a=5.):

    U_p = np.abs(np.diff(u[1:]))
    U_m = np.abs(np.diff(u[:-1]))

    dU_median = np.median(U_p)
    
    i_p = np.where(U_p>dU_median*a)[0]+1
    i_m = np.where(U_p>dU_median*a)[0]+2 
    
    i_nan = np.where(u==np.nan)[0]
    i_error = np.intersect1d(i_p,i_m,assume_unique=True)+i_nan    
    
    c=0.
    for i in i_error:
        c+=1
        u[i]=(u[i-1]+u[i+1])/2    
        
    
    T = 100*c/len(u)
  #  print("percentage of corrected values : "+str(T))
    
    accurate = (T<5.)
    return u,accurate
    
def average(plist):
    keys = plist[0].keys()
    p_mean = {}
    imin = 0
    imax = min([len(p['Xmin']) for p in plist])
   # print(imax)
    
    for key in keys:
        #print(key)
        mat = np.asarray([p[key][imin:imax] for p in plist])
       # print(mat.shape)
        p_mean[key] = np.nanmean(mat,axis=0)
        #print(p_mean[key])
    
    return p_mean
    

def plot(p,tmin,tmax,label='',c=True,fignum=0):
    """
    Plot the position of the vortex as a function time.
    pos is a dictionnary obtained from track.position
    tmin :
        minimum index
    tmax : 
        maximum index
    """
    figs = {}
    keys = ['Xmax','Xmin','Ymin','Ymax']
    subplot = {'Xmax':121,'Xmin':121,'Ymax':122,'Ymin':122}
    
    fig1 = graphes.set_fig(fignum+1)
    fig1.set_size_inches(10,4)
    
    accurate = {key:None for key in keys}
    for key in keys:
      #  print(p[key][tmin:tmax])
        if c:
            p[key][tmin:tmax],accurate[key] = correct(p[key][tmin:tmax],a=5.)
        else:
            accurate[key]=True
            
        if 'Y' in key:
            #print('invert !')
            if np.nanmean(p[key])<0:
                p[key] = -np.asarray(p[key]) #X axis is inverted !
        
        if accurate[key]:   
            graphes.set_fig(fignum+1,subplot[key])
            graphes.graph(p['t'],p[key],fignum=fignum+1,label=label)
            figs.update(graphes.legende('Time (s)',key[0]+' position (mm)',''))
            if 'Y' in key:
                graphes.set_axis(0.05,p['t'][tmax],0,100)
            else:
                graphes.set_axis(0.05,p['t'][tmax],-50,50)
    
    p['d']=np.sqrt((np.asarray(p['Xmin'])-np.asarray(p['Xmax']))**2+(np.asarray(p['Ymin'])-np.asarray(p['Ymax']))**2)
    graphes.graph(p['t'],p['d'][tmin:tmax],fignum=fignum+2,label=label)
    graphes.set_axis(0,p['t'][tmax],0,50)
    figs.update(graphes.legende('Time (s)','Distance (mm)',''))
    
    if accurate['Xmin'] and accurate['Ymin']:
        graphes.graph(p['Ymin'][tmin:tmax],p['Xmin'][tmin:tmax],fignum=fignum+3,label=label)
        figs.update(graphes.legende('X position (mm)','Y position (mm)',''))
        graphes.set_axis(0,60,-50,50)
    
    if accurate['Xmax'] and accurate['Ymax']:
        graphes.graph(p['Ymax'][tmin:tmax],p['Xmax'][tmin:tmax],fignum=fignum+3,label=label)
        
    graphes.graph(p['t'],p['Gammamax'],fignum=fignum+4,label=label)
    figs.update(graphes.legende('Time (s)','Circulation (mm^2/s)',''))
    graphes.set_axis(p['t'][tmin],p['t'][tmax],0,5*10**4)
    
    return figs,accurate    
        
def dispersion(Data,j,savedir):
    tmin = 20
    tmax = 160
    
    accurate = Data['accurate'][j]
    pos = Data['pos'][j]
    Mlist= Data['M'][j]
    A = Data['A'][j]
    
    figs = {}
    graphes
    
    for c,key in enumerate(accurate[0].keys()):
        Y = []
        for i,M in enumerate(Mlist):
            if accurate[i][key]:
                #c+=1
                Y.append(pos[i][key][tmin:tmax])
            else:
                print("not accurate")
        
        Y = np.asarray(Y)
#        print(key)
        t = np.asarray(pos[i]['t'][tmin:tmax])
        graphes.set_fig(1)
        #graphes.errorbar(t,np.nanmean(Y,axis=0),0*t,np.nanstd(Y,axis=0),label='k',fignum=c+4)
#        graphes.graph(t,np.nanstd(Y,axis=0),fignum=1,label=label)
#        graphes.set_axis(0.05,0.45,0,15)
#        figs.update(graphes.legende('Time (s)',key,''))

#        graphes.graph(t**3,np.nanstd(Y,axis=0)**2,fignum=2,label=label)
#        graphes.graph(t**3,1000*t**3,fignum=2,label='r-')
#        graphes.set_axis(0.05**3,0.05,0,70)
#        figs.update(graphes.legende('t^3',key+'^2',''))

        print(A)
        if A==0:
            label='c'
        else:
            label='k-'
            
        if 'X' in key:
            graphes.graphloglog(t,np.nanstd(Y,axis=0),fignum=3,label=label)
            graphes.graphloglog(t,30*t**1.5,fignum=3,label='r--')
          #  graphes.graphloglog(t,8*t**0.5,fignum=3,label='r--')

            graphes.set_axis(0.05,0.6,0.01,50)
            figs.update(graphes.legende('t',key+'',''))

        if 'Y' in key:
#            label='b-'
            graphes.graphloglog(t,np.nanstd(Y,axis=0),fignum=4,label=label)
            graphes.graphloglog(t,30*t**1.5,fignum=4,label='r--')
          #  graphes.graphloglog(t,8*t**0.5,fignum=4,label='r--')
            graphes.set_axis(0.05,0.6,0.01,50)
            figs.update(graphes.legende('t',key+'',''))
            
    #    graphes.graphloglog(t**3,1000*t**3,fignum=2,label='r-')
        #figs.update(graphes.legende('t',key+'',''))
        #print(c)
    return figs

def main(Mlist):
    savedir = title(Mlist[0])    
    
    indices = range(400,2000,10)
    n = len(indices)
    N = len(Mlist)
    print(N)

    figs = {}

    X = np.zeros((N,n,2))
    Y = np.zeros((N,n,2))

    for j,M in enumerate(Mlist):
        print(j)
        Xmin = []
        Xmax = []
        Ymin = []
        Ymax = []
    
        Vmin = []
        Vmax = []
        
        Gammamin = []
        Gammamax = []

        t=[]
        field = 'omega'
        Omega = access.get_all(M,field)
        for i in indices:
            Z = Omega[:,5:,i]
    
            t.append(M.t[i])
    
            xmin,xmax,ymin,ymax = positions(M,i)
        
            vmin,vmax,G = amplitude(M,i)
            Gammamin.append(G[0])
            Gammamax.append(G[1])
  
            Vmin.append(vmin)
            Vmax.append(vmax)
    
            Xmin.append(xmin)
            Xmax.append(xmax)
            Ymin.append(ymin)
            Ymax.append(ymax)
    
        X[j,:,0] = np.asarray(Xmin)
        X[j,:,1] = np.asarray(Xmax)
    
        Y[j,:,0] = np.asarray(Ymin)
        Y[j,:,1] = np.asarray(Ymax)
    
        graphes.graph(t,Xmin,label='bo',fignum=1)
        graphes.graph(t,Xmax,label='ro',fignum=1)
        figs.update(graphes.legende('Time (s)','Horizontal position (mm)',''))
        
        graphes.graph(Xmin,Ymin,label='bo',fignum=2)
        graphes.graph(Xmax,Ymax,label='ro',fignum=2)
        figs.update(graphes.legende('X (mm)','Y(mm)',''))
        
        graphes.graph(t,Vmin,label='b',fignum=3)
        graphes.graph(t,Vmax,label='r',fignum=3)
        figs.update(graphes.legende('Time (s)','Maximum vorticity (s^-1)',''))
        
        graphes.graph(t,Gammamin,label='b',fignum=4)
        graphes.graph(t,Gammamax,label='r',fignum=4)
        figs.update(graphes.legende('Time (s)','Circulation mm^2/s',''))
        graphes.set_axis(0,0.32,-25000,25000)
        
        
    Dx = (X-np.tile(np.nanmean(X,axis=0),(N,1,1)))
    Dy = (Y-np.tile(np.nanmean(Y,axis=0),(N,1,1)))

    D = np.sqrt(Dx**2+Dy**2)

    D_moy = np.nanmean(D,axis=0)
    #print(t)
    #print(D[j,:,1])
#    for j in range(N):
#        graphes.graph(t,D[j,:,0],label='bo',fignum=4)
#        graphes.graph(t,D[j,:,1],label='ro',fignum=4)

    graphes.graph(t,D_moy[:,0],label='b',fignum=6)
    graphes.graph(t,D_moy[:,1],label='r',fignum=6)
    figs.update(graphes.legende('Time (s)','Distance (mm)','Spreading of vortices'))
    graphes.set_axis(0,0.3,0,20)

    graphes.graph(np.power(t,3),np.power(D_moy[:,0],2),label='b',fignum=5)
    graphes.graph(np.power(t,3),np.power(D_moy[:,1],2),label='r',fignum=5)
    graphes.set_axis(0,0.3**3,0,100)

    figs.update(graphes.legende('t^3 (s^3)','d^2 (mm^2)','Spreading of vortices : rescaling'))
    
    
    graphes.save_figs(figs,savedir=savedir,suffix='',prefix='',frmt='pdf',dpi=300,display=True)
    #    fig,axes = panel.make([111],fignum=3+3*j)
    #    fig.set_size_inches(15,15)
    #    i = 1400
    #    graphes.Mplot(M,'enstrophy',i,fignum=3)
    #    graphes.graph(Xmin,Ymin,label='bo',fignum=3+3*j)
    #    graphes.graph(Xmax,Ymax,label='ro',fignum=3+3*j)