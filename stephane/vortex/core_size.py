

import numpy as np
import scipy.ndimage.filters as filters
import scipy.optimize as opt

import os.path

import stephane.display.graphes as graphes
import stephane.analysis.cdata as cdata
import stephane.analysis.vgradient as vgradient
import stephane.jhtd.strain_tensor as strain_tensor

import stephane.manager.access as access

import stephane.tools.browse as browse
import stephane.tools.Smath as Smath

import os
from matplotlib.colors import Normalize
from scipy import ndimage, interpolate, signal
from scipy import optimize as opt
import random
import scipy.io as sio

import stephane.vortex.circulation_core_rad_single_core_stephane as joseph

def compute(M,method='vorticity',fignum=1,display=True):
    """
    Compute the core size measurement using either the 'vorticity' of the 'circulation' method
         vorticity : compute the vorticity field, and fit the resulting distribution by a gaussian bump
         circulation : compute the vorticity, then integrate it over disks of growing size, and
         fit 1-gaussian bump
    INPUT
    -----
    M : Mdata object
        data
    method :
        'vorticity'
        'joseph'
        'circulation'
    OUTPUT
    -----
    """
    if method=='vorticity':
        lc,std_lc = from_vorticity(M,fignum=fignum,display=display)
        print(fignum/2-1,M.dataDir)
        print('Core size : '+str(lc)+' +/- '+str(std_lc)+ ' mm')
        return lc,std_lc,None,None
        
    if method=='circulation':
        lc,std_lc,Gamma,std_Gamma = from_circulation_1(M,fignum=fignum,display=display)
        print('Core size : '+str(lc)+' +/- '+str(std_lc)+ ' mm')
        print('Circulation : '+str(Gamma)+' +/- '+str(std_Gamma)+ ' mm')
        return lc,std_lc,Gamma,std_Gamma
        
    if method=='circulation_2':
        lc,std_lc,Gamma,std_Gamma = from_circulation_2(M,fignum=fignum,display=display)
        print('Core size : '+str(lc)+' +/- '+str(std_lc)+ ' mm')
        print('Circulation : '+str(Gamma)+' +/- '+str(std_Gamma)+ ' mm')
        return lc,std_lc,Gamma,std_Gamma
        
        
    if method=='joseph':
        try:
            lc,std_lc = from_circulation_joseph(M,fignum=fignum,display=display)
        except:
            print('Computation falls down, skip')
            lc,std_lc = 1,1
        print('Core size : '+str(lc)+' +/- '+str(std_lc)+ ' mm')
        return lc,std_lc,None,None
        
    return None

def from_vorticity(M,fignum=1,display=True):        
    ny,nx,nt = M.shape()
    
    field = 'omega'
    M.get(field)
   # M,field = vgradient.compute(M,'vorticity',filter=True)
    Z = getattr(M,field)
    
    x,y = space_axis_vorticity(M)
    
    a=[]
    indices = range(0,nt,1)
    x0 = []
    y0 = []
    index = []
    
    for i in indices:                
        sigma,center = fit_core_size(x,y,Z[...,i],display=False)
        if sigma is not None:
            a.append(sigma)
            index.append(i)
            x0.append(center[0])
            y0.append(center[1])
    
    lc = np.nanmedian(a)    
    std_lc = np.nanstd(a)
    
    #plot example
    error = 0.5
    if len(index)>0 and std_lc/lc<error:
        i_example = index[len(index)//2]
        print('Indice : '+str(i_example))
        sigma,center,figs = fit_core_size(x,y,Z[...,i_example],fignum=fignum,display=True)

        t = np.asarray([M.t[i] for i in index]) 
        graphes.graph(t,a,label='k',fignum=fignum+1)    
#    graphes.graph(t,x0,label='r.',fignum=3)
#    graphes.graph(t,y0,label='b.',fignum=3)
        title = os.path.basename(M.dataDir)
        graphes.set_axis(np.min(t),np.max(t),0.,7.)
        figs.update(graphes.legende('t (s)','Core size',title,display=False))
        save_graphes(M,figs,method='vorticity')
    return lc,std_lc

def from_circulation_1(M,fignum=1,display=True):
    nx,ny,nt = M.shape()
    
    M,field = vgradient.compute(M,'vorticity',filter=True)
    Z = getattr(M,field) #field vorticity : omega
    
    x,y = space_axis_vorticity(M) 
    
    a=[]
    G0=[]
    indices = range(0,nt,1)
    x0 = []
    y0 = []
    index = []

    Omega = getattr(M,'omega')

    start=True
    figs = {}
    for i in indices:
        sigma,A,center = fit_core_circulation(M,fignum=fignum,display=False)              
        #graphes.save_figs(figs,savedir='./Results/'+os.path.basename(M.dataDir)+'/',suffix='_method_circulation_')
        
    #    print(A)
        if sigma is not None:
            a.append(sigma)
            G0.append(A)
            index.append(i)
            x0.append(center[0])
            y0.append(center[1])
            
            if start:
                start=False
                xi = x[center[0]]
                yi = y[center[1]]
                
    r=[np.sqrt((x1-xi)**2+(y1-yi)**2) for x1,y1 in zip(x[x0],y[y0])]
    lc = np.nanmedian(a)    
    std_lc = np.nanstd(a)
    Gamma = np.nanmedian(G0)    
    std_Gamma = np.nanstd(G0)
    
    #plot example
    error = 0.5
    if len(index)>0 and std_lc/lc<error:
        i_example = index[len(index)//2]
#        fit_core_circulation(x,y,Omega[...,i_example],fignum=fignum,display=display)              
        figs.update(graphes.legende('r (mm)','Circulation mm^2/s',''))

        if display:            
            t = np.asarray([M.t[i] for i in index]) 
            graphes.graph(t,a,label='ro',fignum=fignum+2)
            title = os.path.basename(M.dataDir)
            graphes.set_axis(np.min(t),np.max(t),0.,5.)
            figs.update(graphes.legende('t (s)','Core size',title))
        
            graphes.graph(t,G0,label='kp',fignum=fignum+3)
            graphes.set_axis(np.min(t),np.max(t),-20000,0)
            figs.update(graphes.legende('t (s)','Circulation',title))
            
            graphes.graph(t,r,label='b+',fignum=fignum+4)
            graphes.legende('t (s)','Position (mm)',title)
            save_graphes(M,figs,method='circulation')
    
    return lc,std_lc,Gamma,std_Gamma

def from_circulation_2(M,fignum=1,display=True):
    
#    R_list,Gamma,center,factor = compute_circulation_2(M,fignum=fignum)
    lc,G0,center = fit_core_circulation(M,fignum=fignum,display=True)
    
    nx,ny,nt = M.shape()
    
    for i in range(nt):
        R_list,Gamma,center,factor = circulation_2(M,i)
        
        graphes.graph(R_list,Gamma*factor,fignum=fignum,label='k^')
        graphes.legende('r (bmm)','Circulation (mm^2/s)','')
        graphes.set_axis(0,12.,-7000,500)
    
    return None
    
    
def compute_circulation_2(M,fignum=1,display=True):
    #method based on the contour integral (vorticity !)
    nx,ny,nt = M.shape()
        
    for i in range(nt):        
        R_list,Gamma,center,factor = circulation_2(M,i,fignum=fignum)

#        toto
    return R_list,Gamma,center,factor

def circulation_2(M,i,fignum=1,display=False):    
    
    Omega = access.get(M,'omega',i)
    x,y = space_axis_vorticity(M) 
    
    X,Y,data,center,factor = normalize(x,y,Omega[...,0])   

    dx = M.x[0,1]-M.x[0,0]
#print(dx)

    U,d = vgradient.make_Nvec(M,i)   # Z : d+1 dimension np array 

    nx,ny = X.shape
    R_list=np.arange(1.,15.,0.5)
    Gamma = []
    divergence = []
    for b in R_list:
       # print(b)
        tau = strain_tensor.strain_tensor_loc(U,center[0],center[1],d=2,b=b)
        omega,enstrophy = strain_tensor.vorticity(tau,d=2,norm=False)
        div = strain_tensor.divergence_2d(tau,d=2)
        G = (omega[0,0]-div[0,0]) * np.pi*b**2*dx**2
        Gamma.append(G)
        divergence.append(div[0,0]/np.abs(omega[0,0]))

    R_list = np.asarray(R_list) * dx
    
    if display:
        graphes.graph(R_list,Gamma,fignum=fignum,label='bo')
        graphes.legende('r (mm)','Circulation (mm^2/s)','')

        graphes.graph(R_list,divergence,fignum=fignum+1,label='ko')
        graphes.graph(R_list,np.zeros(len(R_list)),fignum=fignum+1,label='r--')
    
        graphes.legende('r (mm)','Relative 2d divergence','')
        graphes.set_axis(0,30*dx,-0.3,0.3)
        
    return R_list,Gamma,center,factor
        

def save_graphes(M,figs,method='vorticity'):
    graphes.save_figs(figs,savedir='./Results/'+os.path.basename(M.dataDir)+'/',suffix='_method_'+method+'_')
    
def compute_circulation(x,y,Z,method='disk'):
    """
    Compute circulation from 
    """
    X,Y,data,center,factor = normalize(x,y,Z)   
    
    # previous method  : integral of vorticity on a disk
    R,theta = Smath.cart2pol(X,Y)
    
    Rsort = np.sort(np.reshape(R,np.prod(R.shape)))
    d = np.diff(Rsort)
    
#    Determine the number of circles to consider for the integral of omega
    dR = np.max(d)
    R_sample = np.arange(0,np.max(R),dR)
#    graphes.graph(Rsort,d,label='k^')   

    dx = np.mean(np.diff(X[0,:]))
    Gamma = np.asarray([np.nansum(data[R<Ri])*dx**2 for Ri in R_sample])
      
    return R_sample,Gamma,center,factor

def space_axis_vorticity(M):
    Z = getattr(M,'omega')
    dimensions = Z.shape
    X = M.x
    Y = M.y
 #   n=0
#    n = (X.shape[0]-dimensions[0])/2
 #   x = X[n:-n,n:-n]
 #   y = Y[n:-n,n:-n]
    return X,Y   
    
def fit_core_size(x,y,Z,fignum=1,display=False):
    """
    Find the half width of a gaussian bump
    INPUT
    -----
    x : 2d np array 
        spatial coordinates (columns)
    y : 2d np array
        spatial coordinates (lines)
    Z : 2d np array
        data to be fitted (typically vorticity field )
    fignum : int
        figure number for the output. Default 1
    display : bool
    OUTPUT
    -----
    a : float
        parameter of the gaussian bump
    center : 2 elements np array
        center coordinates
    """
    ny,nx = Z.shape
    X,Y,data,center,factor = normalize(x,y,Z)   
    R,theta = Smath.cart2pol(X,Y)

    a0=1
    fun = gaussian
    res = opt.minimize(distance_fit,a0,args=(fun,R,data))

    cond = ((center[0]>5) and (center[0]<nx-5) and (center[1]>5) and (center[1]<ny-5))
    if cond:
        a=np.sqrt(res.x)
    else:
        a=None
        
    if display:
        figs = {}
        graphes.graph(R,factor*data,fignum=3,label='ko')
        graphes.graph(R,factor*gaussian(res.x,R),label='r.',fignum=fignum+2)
        graphes.set_axis(0,20,0,factor*1.1)
        figs.update(graphes.legende('r (mm)','Vorticity s^{-1})',''))
        
        graphes.cla(fignum=fignum)
        graphes.color_plot(X,Y,factor*data,fignum=fignum+3)
        graphes.colorbar()
        figs.update(graphes.legende('X (mm)','Y (mm)','Vorticity',display=False,cplot=True))
        return a,center,figs
    else:
        return a,center            
    
#def fit_core_circulation(x,y,Z,fignum=1,display=False):
def fit_core_circulation(M,fignum=1,display=False):
    
    R,Gamma,center,factor = compute_circulation_2(M,fignum=fignum,display=display)
    factor = np.max(Gamma)
    Gamma = Gamma / factor
    #fit_core_circulation(R_list,Gamma,M,fignum=fignum,display=False)
    
    nx,ny,nt = M.shape()
    
    a0=(1,np.max(Gamma))
    fun = circ_gaussian
    res = opt.minimize(distance_fit,a0,args=(fun,R,Gamma))
    
    cond = ((center[0]>5) and (center[0]<nx-5) and (center[1]>5) and (center[1]<ny-5))
    if cond:
        lc=np.sqrt(res.x[0])       
        G0 = factor*res.x[1]
    else:
        lc=None
        G0=None
        
    if display:
        graphes.graph(R,factor*Gamma,fignum=1,label='ko')
        graphes.graph(R,factor*fun(res.x,R),fignum=1,label='r')
        
    return lc,G0,center
    
def normalize(X,Y,Z,sign=False):
    """
    Center and normalize Z around its maximum value (max of |Z|) define on the nodes of grids x and y
    INPUT
    -----
    X : 2d np array
    Y : 2d np array
    Z : 2d np array
    OUTPUT
    -----
    Xshift : 2d np array
        shifted coordinate
    Yshift : 2d np array
        shifted coordinate
    Znorm : 2d np array
        normalized field Z
    center : list with 2 elements
        coordinates of the center in the grid (X,Y)
    factor : normalization factor applied to Z
        can be used to recover the initial amplitude of Z (after fitting operations)
    """
    #large gaussian to find the maximum peak
    Zfilt = filters.gaussian_filter(Z, sigma=5, order=0, output=None)
        
    i = np.argmax([-np.min(Zfilt),np.max(Zfilt)])
    fun = [np.min,np.max]
    
    #Normalize by the maximum amplitude between the min and the max intensity (for signed quantities)
    Zfilt =  Zfilt / fun[i](Zfilt)
    
    factor = fun[i](Z)
    Znorm =  Z / factor
    #print(fun[i].__name__)
        
    i0,j0 = np.unravel_index(np.argmax(Zfilt),Zfilt.shape)
    x0,y0 = (X[0,j0],Y[i0,0])
    
    Xshift = X - x0
    Yshift = Y - y0

#    np.sign(np.nansum(Z))
    
    center = i0,j0
    return Xshift,Yshift,Znorm,center,factor

    
def distance_fit(a,fun,R,Z):
    Z_th = fun(a,R) 
    return np.sum(Z*np.power(Z_th-Z,2))
    
def gaussian(a,R,A=1):
    return A*np.exp(-np.power(R,2)/a)    

def circ_gaussian(a,R):
    return a[1]*(1-gaussian(a[0],R))
    
def from_circulation(M,display=True):
    pass
    
def compile(Mlist,V=None,method='circulation'):
    symbol = {'50':'^','125':'o','250':'s'}
    color = {'circulation':'r','vorticity':'b','joseph':'k'}
    labels = {key:color[method]+symbol[key] for key in symbol.keys()}
    if V==None:
        sub_labels=labels
        piston_v = None
    else:
        piston_v = str(V)
        sub_labels = {piston_v:labels[piston_v]}#,'125':'ro','250':'bs'}
    
    figs={}
    
    for i,M in enumerate(Mlist):        
        piston1 = browse.get_string(M.Sdata.fileCine,'_v',end='.cine',shift=0,display=False,from_end=True)
        piston2 = browse.get_string(M.Sdata.fileCine,'_v',end='_p30mum',shift=0,display=False,from_end=True)
        
        error = 0.25        
        for piston in [piston1,piston2]:
            if piston in sub_labels.keys():
                print(M.dataDir)
                dx = np.mean(np.diff(M.x[0,:]))
                print('Spatial scale : '+str(dx)+' mm/box')
                lc,std_lc,Gamma,std_Gamma=compute(M,method=method,display=False,fignum=(i+1)*2)
                
#                print(piston,dx,lc,std_lc)
                if std_lc/lc<error:
                    graphes.errorbar(dx,lc,[0],std_lc,label=labels[piston],fignum=250)
                    figs.update(graphes.legende('mm/box','Core size (mm)',''))
                    graphes.set_axis(0,1.5,0,6.)
                    
                    if method == 'circulation':
             #   if np.abs(std_Gamma/Gamma)<error:
                        graphes.errorbar(dx,Gamma,[0],std_Gamma,label=labels[piston],fignum=251)
                        figs.update(graphes.legende('mm/box','Circulation (mm^2/s)',''))
                        graphes.set_axis(0,1.5,-2*10**4,0)
                #print(piston,dx,lc,std_lc
                print('')
    
    print('figure',figs)
    print(figs)
    graphes.save_figs(figs,suffix='Compilation_method_'+method+'_v'+piston_v)
   
   
def from_circulation_joseph(M,fignum=1,display=True):
    direc =M.dataDir # NOTE: make sure this is consistent. Some inconsistancy in the file structure/name
    
    num_y,num_x,nt = M.shape()
    
    display_part = False
    
    
    name = direc.split('/')[7] # NOTE: make sure this is consistent. Some inconsistancy in the file structure/name
    frame_begin = 100 # enter the frame range you would like to do the processing on
    frame_end = 120
    frame_by = 1 # How many frames to skip, of form [frame_begin:frame_end:frame_by]
    ppmm = 1/M.fx # pixels per mm
    radii = np.arange(10,250,15) # radii to calculate the circulation at (in pixels)
    num_pts = 800 # number of points to have on the circles
    position = 'in' # which side of the core is being looked at. in = negative vorticity, out = positive

    # all automatic
    file_dirs = sorted(glob.glob(os.path.join(direc,'*.txt'))) # finding the names of the files data is kept in
    fs = float(name.split('_')[6][3:])
    result_dir = './'+name+'_circ_to_core_DPIV/' # name of the result directory.
    if not os.path.isdir(result_dir):
        os.mkdir(result_dir)
    core_radii = {} # Dictionary which will equate a frame number to the core radius at that frame
    rainbow_cm = cm.gist_rainbow
    rad_list = [] # list of radii found
    if not os.path.isdir(os.path.join(result_dir,'%s_figs/'%position)): # Directory to save the figures in. Is subdir of the result directory
        os.mkdir(os.path.join(result_dir,'%s_figs/'%position))

    # Get a first estimate of the core radius to make a gaussian bump to convolve with the other images
    # Picks out three random frames which are used to get an estimate of the core size.
    rad_est = []
    options = [random.randrange(frame_begin,frame_end+1),random.randrange(frame_begin,frame_end+1),random.randrange(frame_begin,frame_end+1)] # 3 random frames
    for frame_no in options: # need to slice these so we only take frames where it is actually in the frame

        # unload all of the information
        x_array, y_array, ux_array, uy_array = np.loadtxt(file_dirs[frame_no],skiprows=3,usecols=(0,1,2,3),unpack=True,delimiter=',')
        frame = frame_no
        x_array = x_array.astype(int)
        y_array = y_array.astype(int)

        # finding information about which points were taken in PIV
        x_beg = x_array[0]
        x_end = x_array[-1]
        y_beg = y_array[0]
        y_end = y_array[-1]
        y_by = y_array[1] - y_array[0]
        num_y = (y_end - y_beg)/y_by + 1
        x_by = (x_end - x_beg)/(ux_array.shape[0]/(num_y)-1)
        num_x =(x_end - x_beg)/x_by + 1
        
         # NOTE: This is particular. I'm not sure of an easy way to find this when all of them are just in a row. Would probably have to do length of differneces to find a number then multiply
       # print(x_by)
        piv_roi = [0,x_end,0,y_end]
        # Make the arrays 2D        
        plane_y, plane_x = np.mgrid[y_beg:y_end+1:y_by,x_beg:x_end+1:x_by]
        plane_ux = M.Ux[...,frame_no]#np.transpose(np.reshape(ux_array,(num_x,num_y)))
        plane_uy = M.Uy[...,frame_no]#np.transpose(np.reshape(uy_array,(num_x,num_y))) # NOTE: x goes from right to left!!!
        plane_ux.shape

        # removing NaN and spurious points from the data
        good_plane_ux = joseph.rm_nans(plane_ux)
        good_plane_uy = joseph.rm_nans(plane_uy)
        best_plane_ux, best_plane_uy = joseph.rm_big(good_plane_ux,good_plane_uy)
        best_ux_array = np.transpose(best_plane_ux).flatten()
        best_uy_array = np.transpose(best_plane_uy).flatten()

        # compute curl
        duy_dx = joseph.partial_deriv(best_plane_uy,'x',x_by) # NOTE Changed this so no -1 in the dx direction
        dux_dy = joseph.partial_deriv(best_plane_ux,'y',y_by)
        plane_curl = duy_dx - dux_dy

        if position == 'out':
            # take maximum of vorticity as center of vortex core
            ext_curl = np.where(plane_curl == plane_curl.max())
            ext_cx = ext_curl[1]
            ext_cy = ext_curl[0]
        else:
            # take minimum of vorticity as center of vortex core
            ext_curl = np.where(plane_curl == plane_curl.min())
            ext_cx = ext_curl[1] # NOTE these are wrt the PIV index NOT the real image pixels
            ext_cy = ext_curl[0]

        # interpolate
        i_ux = joseph.interp_flow_component(x_array,y_array,best_ux_array,piv_roi) # returning a grid interpolated to every 1 (real) pixel instead of every 16 which is what we get from PIV
        i_uy = joseph.interp_flow_component(x_array,y_array,best_uy_array,piv_roi)
        phis = np.arange(0,2*np.pi,2*np.pi/float(num_pts)) # for plotting circles
        gammas = []
        for rad in radii:
            circx, circy = joseph.circPts(rad,phis,ext_cx*x_by,ext_cy*y_by) # circular path to do integral around
            ds = rad*2*np.pi/float(num_pts) # element of arc length
            i_ux_path = ndimage.map_coordinates(i_ux,[circy,circx],order=1) # points of ux interpolated on the path
            i_uy_path = ndimage.map_coordinates(i_uy,[circy,circx],order=1)
            gamma = sum(i_uy_path * np.cos(phis) - i_ux_path * np.sin(phis))*ds # the integral
            gammas.append(gamma)
        gammas = np.asarray(gammas) # found circulations
        real_gammas = gammas * fs * (1/ppmm)**2 # converting to real units, mm^2/s
        real_radii = radii / ppmm # converting to real units, mm
        popt, pcov = opt.curve_fit(joseph.circ_by_rad,real_radii,real_gammas,p0=(real_gammas.max(),4.)) # fitting to the correct distribution, give rough first est of core size of 4 mm
        sigma = np.abs(popt[1]) # should be core size in mm
        rad_est.append(sigma) # just need an estimate
  #      print('done with estimate %d'%frame_no)
    med_rad_est = np.median(rad_est)*ppmm # estimate for mean radius from the 3 random frames
    print('estimated radius: %.2f mm'%(med_rad_est/ppmm))

    # need to make the thing to convolve
    # making the thing to convolve
    y,x = np.mgrid[0:num_y,0:num_x]
    gauss_bump = joseph.ellip_gauss_func(x,y,num_x/2.,num_y/2.,med_rad_est/x_by,med_rad_est/y_by)

    # this is the real processing now
    k = frame_begin
    for file_dir in file_dirs[frame_begin:frame_end+1:frame_by]: # need to slice these so we only take frames where it is actually in the frame

        # unload all of the information
        x_array, y_array, ux_array, uy_array = np.loadtxt(file_dir,skiprows=3,usecols=(0,1,2,3),unpack=True,delimiter=',')
        frame = k
        x_array = x_array.astype(int)
        y_array = y_array.astype(int)

        # finding information about which points were taken in PIV
        x_beg = x_array[0]
        x_end = x_array[-1]

        y_beg = y_array[0]
        y_end = y_array[-1]
        y_by = y_array[1] - y_array[0]
        num_x =(x_end - x_beg)/x_by + 1
        num_y = (y_end - y_beg)/y_by + 1
        piv_roi = [0,x_end,0,y_end]

        # Make the arrays 2D
        plane_y, plane_x = np.mgrid[y_beg:y_end+1:y_by,x_beg:x_end+1:x_by]
        plane_ux = np.transpose(np.reshape(ux_array,(num_x,num_y)))
        plane_uy = np.transpose(np.reshape(uy_array,(num_x,num_y))) # NOTE: x goes from right to left!!!

        # removing NaN and spurious points from the data
        good_plane_ux = joseph.rm_nans(plane_ux)
        good_plane_uy = joseph.rm_nans(plane_uy)
        best_plane_ux, best_plane_uy = joseph.rm_big(good_plane_ux,good_plane_uy)
        best_ux_array = np.transpose(best_plane_ux).flatten()
        best_uy_array = np.transpose(best_plane_uy).flatten()

        # compute curl
        duy_dx = joseph.partial_deriv(best_plane_uy,'x',x_by) # NOTE Changed this so no -1 in the dx direction
        dux_dy = joseph.partial_deriv(best_plane_ux,'y',y_by)
        plane_curl = duy_dx - dux_dy

        # convolve with a gaussian bump to find vortex core
        convolution = signal.convolve2d(gauss_bump,plane_curl)

        # find and rescale the convolution center
        if position == 'out':
            conv_ext_indx = np.where(convolution == convolution.max()) # IDK Why but seemed like we had to go to -1 instead of +1 here
            conv_ext_x = (conv_ext_indx[1][0] - num_x-1)*x_by + x_beg
            conv_ext_y = (conv_ext_indx[0][0] - num_y-1)*y_by + y_beg
            # convolution only gives us the positions relative to the center of the image so translate
            ext_cx = conv_ext_x + x_end/2.
            ext_cy = conv_ext_y + y_end/2.
        else:
            conv_ext_indx = np.where(convolution == convolution.min())
            conv_ext_x = (conv_ext_indx[1][0] - num_x-1)*x_by + x_beg
            conv_ext_y = (conv_ext_indx[0][0] - num_y-1)*y_by + y_beg
            # convolution only gives us the positions relative to the center of the image so translate
            ext_cx = conv_ext_x + x_end/2.
            ext_cy = conv_ext_y + y_end/2.

        # sub-pixel interpolation for convolution
        ygrid, xgrid = np.mgrid[0:num_y,0:num_x]
        xgrid_f = xgrid.ravel()*x_by # should be in real pixel location now
        ygrid_f = ygrid.ravel()*y_by
        if position == 'out':
            init_guess = (ext_cx,ext_cy, med_rad_est,med_rad_est,1)
        else:
            init_guess = (ext_cx,ext_cy, med_rad_est,med_rad_est,-1)
        popt, pcov = opt.curve_fit(joseph.gauss_func_cfit, (xgrid_f,ygrid_f), plane_curl.ravel(),p0 = init_guess)
        ext_cx_n, ext_cy_n = popt[:2] # should be an even better approximation

        # interpolate
        i_ux = joseph.interp_flow_component(x_array,y_array,best_ux_array,piv_roi) # returning a grid interpolated to every 1 (real) pixel instead of every ?? which is what we get from PIV
        i_uy = joseph.interp_flow_component(x_array,y_array,best_uy_array,piv_roi)
        phis = np.arange(0,2*np.pi,2*np.pi/float(num_pts)) # for circles
       # print 'frame %d'%frame

        pos_dir = '%s_figs/'%position
        gammas = []
        for rad in radii:
            circx, circy = joseph.circPts(rad,phis,ext_cx_n,ext_cy_n) # circular path to do integral around
            ds = rad*2*np.pi/float(num_pts) # element of arc length
            i_ux_path = ndimage.map_coordinates(i_ux,[circy,circx],order=1) # points of ux interpolated on the circular path
            i_uy_path = ndimage.map_coordinates(i_uy,[circy,circx],order=1)
            gamma = sum(i_uy_path * np.cos(phis) - i_ux_path * np.sin(phis))*ds # the integral
            gammas.append(gamma)
        gammas = np.asarray(gammas) # found circulations
        real_gammas = gammas * fs * (1/ppmm)**2 # converting to real units, mm^2/s
        real_radii = radii / ppmm # converting to real units, mm
        popt, pcov = opt.curve_fit(joseph.circ_by_rad,real_radii,real_gammas,p0=(real_gammas.max(),med_rad_est)) # fitting to the correct distribution
        
        if display_part:
            graphes.set_figure(1)
            fig = plt.figure(figsize=(18.0,8.0)) # makes the figure 18x8 inches~
            plt.subplot(221)
            plt.title('ux for frame %d'%frame)
            norm = joseph.MidpointNormalize(midpoint=0) # for making sure white is really 0 in the color map
            plt.imshow(best_plane_ux,norm=norm,cmap=cm.bwr,extent=[0,x_end,0,y_end]) # plot of ux

        # plotting circles we took the integrals over
            for rad in radii:
                circx, circy = joseph.circPts(rad,phis,ext_cx_n,ext_cy_n)
                circy = y_end - circy # because the origin of the image is at the top but origin for plot is at the bottom
                plt.plot(circx,circy,color=rainbow_cm(float(rad)/radii.max()),zorder=10)
            plt.subplot(222)
            plt.title('uy for frame %d'%frame)
            norm = joseph.MidpointNormalize(midpoint=0) # for making sure white is really 0 in the color map
            plt.imshow(best_plane_uy,norm=norm,cmap=cm.bwr,extent=[0,x_end,0,y_end]) # plot of uy

        # plotting circles we took the integrals over
        for rad in radii:
            circx, circy = joseph.circPts(rad,phis,ext_cx_n,ext_cy_n)
            circy = y_end - circy # because the origin of the image is at the top but origin for plot is at the bottom
            plt.plot(circx,circy,color=rainbow_cm(float(rad)/radii.max()))

        # plotting calculated circulation as a function of radius of path
        if display_part:
            plt.subplot(212)
            fit_dom = np.linspace(0,real_radii.max(),100)
            fit_range = joseph.circ_by_rad(fit_dom,popt[0],popt[1])
            plt.plot(fit_dom,fit_range,color='k',zorder=1)
            for i, rad in enumerate(real_radii):
                plt.scatter(rad,real_gammas[i],marker='o',color=rainbow_cm(float(rad)/real_radii.max()),zorder=2)
            plt.xlabel('radius (mm)')
            plt.ylabel('circulation (mm^2/s)')
            plt.xlim(0,real_radii.max()*1.01)
            plt.ylim(real_gammas.min()-np.abs(real_gammas.min())*0.01,real_gammas.max()+np.abs(real_gammas.max())*0.01)
            plt.savefig(os.path.join(result_dir,pos_dir,'circulation_v_rad_frame_%d.png'%frame))
            plt.close()
        
        sigma = np.abs(popt[1]) # should be core size in mm
        core_radii[frame] = sigma # saving    
        rad_list.append(sigma)
        k += frame_by

#        print 'measured core radius (%s): %.2f'%(position,sigma)

        # plotting vorticity with an overlay of the core
        if display_part:
            norm = joseph.MidpointNormalize(midpoint=0) # for making sure white is really 0 in the color map
            plt.imshow(plane_curl,norm=norm,cmap=cm.bwr,extent=[0,x_end,0,y_end])
            circx, circy = joseph.circPts(sigma*ppmm,phis,ext_cx_n,ext_cy_n)
            circy = y_end - circy # because the origin of the image is at the top but origin for plot is at the bottom
            plt.plot(circx,circy,'g')
            plt.scatter([ext_cx_n],[y_end - ext_cy_n],marker='x',color='g')
            plt.title('frame %d'%frame)
            plt.savefig(os.path.join(result_dir,pos_dir,'found_core_overlay_frame_%d.png'%frame))
            plt.close()
            

    # plotting core radius over time and the mean
    t = np.arange(0,len(rad_list))*frame_by
    mean_rad = sum(rad_list)/len(rad_list)
    median_rad = np.median(rad_list)
    top_y = max(rad_list)
    
    if display_part:
        fig = plt.figure(figsize=(18.0,8.0)) # makes the figure 18x8 inches~ (basically able to save it maximized now!)
        plt.plot(t,rad_list,'b.-')
        plt.plot([t[0],t[-1]],[median_rad,median_rad],'r')
        plt.title('%s measured core radius over time'%position)
        plt.xlabel('frame number from beginning of analysis')
        plt.ylabel('core radius (mm)')
        plt.ylim(0,max(rad_list)*1.01)
        plt.xlim(t[0]*0.99,t[-1]*1.01)
        plt.savefig(os.path.join(result_dir,'core_rad_time.png'))
        plt.close()

    # saving the info we collected
    core_rad_file = open(os.path.join(result_dir,'core_radii.p'),'wb')
    cPickle.dump(core_radii,core_rad_file)
    core_rad_file.close()

    return median_rad,np.std(rad_list)
