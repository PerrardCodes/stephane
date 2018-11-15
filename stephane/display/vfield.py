# -*- coding: utf-8 -*-
"""
Created on Wed Oct 21 15:39:23 2015

@author: stephane
"""

import numpy as np
import stephane.display.graphes as graphes
import os.path
import matplotlib.pyplot as plt

import stephane.jhtd.strain_tensor as strain_tensor
import stephane.analysis.vgradient as vgradient

def make_2dmovie(M,name,Range=None,fignum=1,local=True,log=False,vmin=0,vmax=1,filt=False):
    """
    Movie of the colormap of the velocity modulus U over time  
    INPUT
    -----
    M : Mdata class instance, or any other object that contains the following fields :
        methods : shape()
        attributes : Ux, Uy
        Sdata object
        Ids object
        ...
    name : name of the field to be plotted. example : Ux,Uy, vorticity, strain
    Range : np array
    fignum : int. default value is 1
    Dirname : string
        Directory to save the figures
    log : bool. default value is True
    OUTPUT
    -----
    None
    """
    if Range is not None:
        start,stop,step = tuple(Range)
    else:
        start,stop,step = tuple([0,M.shape()[-1],5])
#    U=np.sqrt(S.Ux**2+S.Uy**2)
    #by default, pictures are save here :
    if local:
        Dirlocal = '/Users/stephane/Documents/Experiences_local/Accelerated_grid'
    else:
        Dirlocal = os.path.dirname(M.Sdata.fileCine)
    Dirname = name
    if filt:
        Dirname = Dirname + '_filtered'
    
    Dir =  Dirlocal + '/PIV_data/'+M.Id.get_id()+'/'+Dirname+'/'
    print(Dir)
    
    fig = graphes.set_fig(fignum)
    graphes.set_fig(0) #clear current figure    
    
    print('Compute '+name)
    M,field = vgradient.compute(M,name,filter=filt)

    for i in range(start,stop,step):
        #Z = energy(M,i)
        graphes.Mplot(M,field,i,fignum=fignum,vmin=vmin,vmax=vmax,log=log,auto_axis=True)
#        graphes.color_plot(Xp,Yp,-Z,fignum=fignum,vmax=700,vmin=0)
        
        if i==start:
            cbar = graphes.colorbar()#fignum=fignum,label=name+'(mm^2 s^{-2})')        
        else:
            print('update')
           # cbar.update_normal(fig)

        filename = Dir + 'V' + str(i)
        graphes.save_fig(fignum,filename,frmt='png')
        
        
def energy(M,i=None):
    """
    Return the kinetic energy of 2d/3d piv field
    INPUT
    -----
    M : Mdata object
        with attributes Ux, Uy (d+1 np array)
    i : int. index of time axis (last axis)
    """
    if i is not None:
        E=M.Ux[...,i]**2+M.Uy[...,i]**2
    else:
        E=M.Ux**2+M.Uy**2
    return E
        
def compute_vorticity(M,i):
    """
    
    """
    U = M.Ux[...,i]
    V = M.Uy[...,i]
    
    X = M.x
    Y = M.y

    x = M.x[0,:]
    y = M.y[:,0]
    
    n=3
    Xp,Yp = np.meshgrid(x[n:-n],y[n:-n])
    
    dx=np.mean(np.diff(x))
    
    dimensions = np.shape(U)+(2,)
    Z = np.reshape(np.transpose([V,U],(1,2,0)),dimensions)
    dZ = strain_tensor.strain_tensor(Z,d=2)
    omega,enstrophy = strain_tensor.vorticity(dZ,d=2,norm=False)
    
    omega = omega/dx
    
    return Xp,Yp,omega    
    
def plot(x,y,U,fignum=1,vectorScale=10**8):
    """
    Plot a 2d velocity fields with color coded vectors
    Requires fields for the object M : Ux and Uy
    INPUT
    -----	
    M : Mdata set of measure 
    frame : number of the frame to be analyzed
    fignum (opt) : asking for where the figure should be plotted
    
    OUTPUT
    ------
    None
    	"""
    Ux=U[:,0]
    Uy=U[:,1]
    
    colorCodeVectors=False
    refVector = 1.
#    vectorScale = 100
    vectorColormap='jet'
    
    #bounds
    #chose bounds from the histograme of E values ?
    scalarMinValue=0
    scalarMaxValue=100
    
    # make the right figure active
    graphes.set_fig(fignum)
    
    # get axis handle
    ax = plt.gca()
    ax.set_yticks([])
    ax.set_xticks([])

    E = np.sqrt(Ux**2 + Uy**2)
    Emoy = np.nanmean(E)
    
    if colorCodeVectors:
        Q = ax.quiver(x,y,Ux/Emoy,Uy/Emoy, E, \
                scale=vectorScale/refVector, 
                scale_units='width', 
                cmap=plt.get_cmap(vectorColormap), 
                clim=(scalarMinValue, scalarMaxValue),
                edgecolors=('none'),
                zorder=4)
    #elif settings.vectorColorValidation:
    #    v = 1
    #    #ax.quiver(x[v==0], y[v==0], ux[v==0], uy[v==0], \
        #    scale=vectorScale/refVector, scale_units='width', color=[0, 1, 0],zorder=4)
    #    Q = ax.quiver(x[v==1], y[v==1], ux[v==1], uy[v==1], \
    #                  scale=vectorScale/refVector, scale_units='width', color='red',zorder=4)
    else:
        Q = ax.quiver(x,y,Ux/E,Uy/E, scale=vectorScale/refVector, scale_units='width', zorder=4) #, color=settings.vectorColor
    
    graphes.legende('$x$ (mm)','$y$ (mm)','')
    
        # add reference vector
    #if settings.showReferenceVector:
    #        plt.quiverkey(Q, 0.05, 1.05, refVector, str(refVector) + ' m/s', color=settings.vectorColor)

    #overwrite existing colorplot
    graphes.refresh(False)

####### to add :

#   1. streamline plots
#   2. Strain maps
#   3. Vorticity maps

#for i in range(10,5000,1):
#    vfield_plot(M_log[4],i,1)
   # input()
   
#    strain = np.sqrt(np.power(dZ[...,0,0],2)+np.power(dZ[...,1,1],2))
#    graphes.color_plot(Xp,Yp,strain,fignum=7)
#    graphes.colorbar()
#    graphes.legende('X','Y','Strain')
    
    
#    graphes.color_plot(X,Y,R,fignum=3,vmin=0,vmax=10)
    
#    print(omega)

    
    