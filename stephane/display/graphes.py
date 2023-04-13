

import matplotlib.pyplot as plt
import numpy as np
import os
import pickle

from math import *
#import stephane.tools.fitting as fitting
import stephane.hdf5.h5py_s as h5py_s

def plot(fun,x,y,fignum=1,label='-',subplot=None,**kwargs):
    """
    plot a graph using the function fun
    fignum can be specified
    any kwargs from plot can be passed
    Use the homemade function refresh() to draw and plot the figure, no matter the way python is called (terminal, script, notebook)
    """
    #set_fig(fignum,subplot=subplot)
    fun(x,y,label,**kwargs)
    refresh()   

def graph(x,y,fignum=1,label='-',subplot=None,**kwargs):
    """
    plot a graph using matplotlib.pyplot.plot function
    fignum can be specified
    cut x data if longer than y data
    any kwargs from plot can be passed
    Use the homemade function refresh() to draw and plot the figure, no matter the way python is called (terminal, script, notebook)
    """    
    xp = np.asarray(x)
    yp = np.asarray(y)
    if len(xp)>len(yp):
        print("Warning : x and y data do not have the same length")
        xp=xp[:len(yp)]    
        
    plot(plt.plot,xp,yp,fignum=fignum,label=label,subplot=subplot,**kwargs)
    
    
def graphloglog(*args,**kwargs):
    plot(plt.loglog,*args,**kwargs)

def semilogx(*args,**kwargs):
    plot(plt.semilogx,*args,**kwargs)

def semilogy(*args,**kwargs):
    plot(plt.semilogy,*args,**kwargs)
    
def errorbar(x,y,xerr,yerr,fignum=1,label='k^',subplot=None,**kwargs):
    """
    plot a graph using matplotlib.pyplot.errorbar function
    fignum can be specified
    cut x data if longer than y data
    any kwargs from plot can be passed
    """
    set_fig(fignum,subplot=subplot)        
    plt.errorbar(x,y,yerr,xerr,label,**kwargs)
    refresh()

def set_fig(fignum,subplot=None):
    if fignum==-2:
        fig=None
        pass
    if fignum==-1:
        fig=plt.figure()
    if fignum==0:
        fig=plt.cla()
    if fignum>0:
        fig=plt.figure(fignum)
        
    if subplot is not None:
        # a triplet is expected !
        ax = fig.add_subplot(subplot)
        return fig,ax
    return fig
    
def cla(fignum):
    set_fig(fignum)
    plt.cla()
    
def set_axis(xmin,xmax,ymin,ymax):
    #set axis
    plt.xlim(xmin=xmin,xmax=xmax)
    plt.ylim(ymin=ymin,ymax=ymax)        
    
def time_label(M,frame):
    Dt = M.t[frame+1]-M.t[frame]
    title = 't = '+str(floor(M.t[frame]*1000)/1000.)+' s, Dt = '+str(floor(Dt*10**4)/10.)+' ms'
    return title
    
def pause(time=3):
    plt.pause(time)
    
def refresh(hold=True,block=False,ipython=True):
    """
    Refresh the display of the current figure. 

    INPUT 
    hold (opt): False if the display has overwritten.

    OUTPUT 
    None
    """    
    if not ipython:           
#        plt.pause(0.1)
#        plt.draw()   
        plt.hold(hold)
        plt.show(block)

def subplot(i,j,k):
    plt.subplot(i,j,k)

def make_title(keylist,d):
    s = r''
    for key in keylist:
        val = str(np.round(d[key]['value'],decimals=1))
        s = s+r'$'+d[key]['legend']+'$ = '+val+d[key]['unit']+', '
    s = s[:-2]
    return s

def legende(x_legend,y_legend,title,ax=None,font=18,display=False,cplot=False,show=True,**kwargs):
    """
    Add a legend to the current figure
        Contains standard used font and sizes
        return a default name for the figure, based on x and y legends
    INPUT
    -----
    x_legend : str
        x label
    y_legend : str
        y label
    title : str
        title label
    colorplot : bool
        default False
        if True, use the title for the output legend
    OUTPUT
    -----
    fig : dict
        one element dictionnary with key the current figure number
        contain a default filename generated from the labels
    """
    if ax is None:
        ax = plt.gca()
    else:
        plt.sca(ax)
    #additionnal options ?
    plt.rc('font',family='Times New Roman')
    plt.rc('text', usetex=False)
    
    if r'\frac' in x_legend:
        plt.xlabel(x_legend,fontsize=font+4)
    else:
        plt.xlabel(x_legend,fontsize=font)
        
    if r'\frac' in y_legend:
        plt.ylabel(y_legend,fontsize=font+4,rotation = 0)
        ax.yaxis.set_label_coords(-0.13, 0.5) 
    else:
        plt.ylabel(y_legend,fontsize=font,*kwargs)#,rotation = 0)        
#    plt.ylabel(y_legend,fontsize=font)
    plt.title(title,fontsize=font)
    
    if show:
        refresh()

    #rec
    fig = figure_label(x_legend,y_legend,title,display=display,cplot=cplot) #fig is a dictionnary where the key correspond to the fig number and the element to the automatic title
    fig = get_data(fig)

    return fig

def get_data(fig,cplot=False):
    current = plt.gcf()
    lines = plt.gca().get_lines()
    
    Dict = {}
    for i,line in enumerate(lines):
        xd = line.get_xdata()
        yd = line.get_ydata()
        Dict['xdata_'+str(i)]=xd
        Dict['ydata_'+str(i)]=yd
        
        if cplot:
            zd = line.get_zdata()
            Dict['zdata'+str(i)]=zd
                
    fig[current.number]['data']=Dict  
    return fig
    
def figure_label(x_legend,y_legend,title,display=True,cplot=False,include_title=False):    
    #generate a standard name based on x and y legend, to be used by default as a file name output
    x_legend = remove_special_chars(x_legend)
    y_legend = remove_special_chars(y_legend)
    title = remove_special_chars(title)
    
    fig={}
    current = plt.gcf()
    fig[current.number]={}
    if cplot:
        fig[current.number]['fignum']=title#+'_'+x_legend+'_'+y_legend #start from the plotted variable (y axis)
    else:
        fig[current.number]['fignum']=y_legend+'_vs_'+x_legend #start from the plotted variable (y axis)
        
    if include_title:
        fig[current.number]['fignum']=y_legend+'_vs_'+x_legend+'_'+title
        
    if display:
        print(current.number,fig[current.number])
    return fig
 
def remove_special_chars(string,chars_rm=['$','\ ','[',']','^','/',') ','} ',' '],chars_rp=['{','(',',','=','.']):
    """
    Remove characters from a typical latex format to match non special character standards
    INPUT
    -----
    string : str
        input string
    chars_rm : str list. Default value : ['$','\ ',') ']
        char list to be simply removed
    chars_rp : str list. Default value : ['( ',',']
        char list to be replaced by a '_'
    OUTPUT
    -----
    string : str
        modified string
    """
    for char in chars_rm:
        string = string.replace(char[0],'')
    for char in chars_rp:
        string = string.replace(char[0],'_')
    return string
    
def title(M):
    #standard title format to know from what date is has been created
    #need something much more sophisticated than that !
    #create a dictionnary of keyword to be added in the title, and add every element into a string formatted style
    tdict={}
    tdict['type']=M.param.typeplane
    if M.param.typeplane=='sv':
        tdict['X']=M.param.Xplane
    if M.param.typeplane=='bv':
        tdict['Z']=M.param.Zplane
    
    tdict['fx']=M.param.fx
    
    title = M.Id.get_id()+','
    for key in tdict.keys():
        title += key + '=' + str(tdict[key])+','

    return title    
    
def colorbar(fignum=-2,label=''):
    set_fig(fignum)
    c=plt.colorbar()    
    
    c.set_label(label,fontsize=18)
    return c
    
def set_title(M,opt=''):
    #if Zplane attribute exist !
    title='Z= '+str(int(M.param.Zplane))+', '+M.param.typeview+', mm, '+M.Id.get_id()+', '+opt 
    plt.title(title,fontsize=18)
    return title

def clegende(c,c_legend):
    c.set_label(c_legend)    

def save_graphes(M,figs,prefix='',suffix=''):
    save_figs(figs,savedir='./Results/'+os.path.basename(M.dataDir)+'/',prefix=prefix,suffix=suffix)

def save_figs(figs,savedir='./',suffix='',prefix='',frmt='pdf',dpi=300,display=False):
    """
    save a dictionnary of labeled figures using dictionnary elements 
        dict can be autogenerated from the output of the graphes.legende() function
        by default the figures are saved whithin the same folder from where the python code has been called
    INPUT
    -----
    figs : dict of shape {int:str,...}
        the keys correspond to the figure numbers, the associated field
    savedir : str. default : './'
    frmt : str
        file Format
    dpi : int
        division per inch. set the quality of the output
    OUTPUT
    -----
    None
    """
    c=0
    filename=''
    for key in figs.keys():
        fig = figs[key]        
        #save the figure
        filename = savedir+prefix+fig['fignum']+suffix
        save_fig(key,filename,frmt=frmt,dpi=dpi)
        c+=1

        #save the data
        if not fig['data']=={}:
            with open(filename+'.pkl', 'wb') as handle:
                pickle.dump(fig['data'], handle, protocol=pickle.HIGHEST_PROTOCOL)
            #print(fig['data'])
            #try:
            #    h5py_s.save(filename,fig['data'])
            #except:
            #    print('Data not saved')
    if display:
        print('Number of auto-saved graphs : '+str(c))
        print(filename)
    
def save_fig(fignum,filename,frmt='pdf',dpi=300,overwrite=False):
    """
    Save the figure fignumber in the given filename
    INPUT
    -----
    fignumber : number of the fig to save
    filename : name of the saving file
    fileFormat (opt)
    dpi (opt) : number of dpi for image based format. Default is 300
    OUTPUT 
    -----
    None
    """
    if os.path.dirname(filename)[0]=='.':
        Dir = '/Users/stephane/Documents/git/stephane/stephane/iPython_notebooks'+os.path.dirname(filename)[1:]
        filename = Dir +'/'+ os.path.basename(filename)
    else:
        Dir = os.path.dirname(filename)
    if not os.path.isdir(Dir):
        os.makedirs(Dir)

    filename=filename+'.'+frmt
    if fignum!=0:
        plt.figure(fignum)    
    
    if not os.path.isfile(filename):
        plt.savefig(filename, dpi=dpi)
    elif overwrite:
        plt.savefig(filename, dpi=dpi)
    else:
        print("figure already exists")
        
def plot_axes(fig,num):
    ax = fig.add_subplot(num)
    ax.set_aspect('equal', adjustable='box')
    #graphes.legende('','','Front view')
    #draw_fieldofview(M.Sdata,ax3,view='front')
    return ax
            
    
def color_plot(x,y,Z,fignum=1,vmin=0,vmax=0,log=False,show=False,cbar=False):
    """
    Color coded plot
    INPUT
    -----	
    x : 2d numpy array
    y : 2d numpy array
    Z : 2d numpy array 
    OUTPUT
    ------
    None
    	"""
    fig,ax = set_fig(fignum,subplot=111)
    #ax.axis('off')
        
    if log:
        Z = np.log10(Z)
         
   # vmin=-3
 #   vmax=max([1,log(np.std(Z)*2)])
  #  Z=np.log(Z)
#    ax = plot_axes(fig,111)
  #  ax.axis('off')
    if vmin==vmax==0:
        c=plt.pcolormesh(x,y,Z)#,shading='gouraud')
    else:
        c=plt.pcolormesh(x,y,Z,vmin=vmin,vmax=vmax)
    #set_axis(np.min(x),np.max(x),np.min(y),np.max(y))
#    plt.axis((np.min(x),np.max(x),np.min(y),np.max(y)))
#    plt.xlim(xmin=np.min(x),xmax=np.max(x))
#    plt.ylim(xmin=np.min(x),xmax=np.max(x))
    if cbar:
        colorbar()
    if show:
        refresh()
  #  return fig,ax,c
    
def get_axis_coord(M,direction='v'):
    X = M.x
    Y = M.y
    
    if hasattr(M,'param'):
        if M.param.angle==90:
            Xb = X
            Yb = Y
            X = Yb
            Y = Xb
        
    return X,Y    
#    return rotate(X,Y,M.param.angle)
    
def rotate(X,Y,angle):
    angle = angle/180*np.pi
    return X*np.cos(angle)-Y*np.sin(angle),Y*np.cos(angle)+X*np.sin(angle)
    
def Mplot(M,field,frame,auto_axis=False,step=1,W=None,Dt=None,fignum=1,show=False,vmin=0,vmax=0,log=False,display=False,tstamp=False,compute=False,cbar=False,colorbar=False):

    import stephane.pprocess.check_piv as check
    import stephane.manager.access as access

    data = access.get(M,field,frame,step=1,compute=compute)
    dimensions = data.shape
    
    if field=='strain':
        #tensor variable. chose the trace (2d divergence !):
        data = data[...,1,1,:]+data[...,0,0,:]
        #print(data)
        
    X,Y = get_axis_coord(M)
    jmin = 0
    data = data[:,jmin:]
    X = X[:,jmin:]
    Y = Y[:,jmin:]
    
    t = M.t[frame]          
    ft = M.t[frame+1]-M.t[frame]
    dx = np.mean(np.diff(M.x[0,:]))
    
    if dx==0:
        dx=1
    
    if vmin==0 and vmax==0:
        if auto_axis:
            std = np.sqrt(np.nanmedian(np.power(data,2)))
            vmax = 10*std
            vmin = -vmax
            if field in ['E','enstrophy']:
                vmin=0

            else:
                if W is None:
                    vmin,vmax = check.bounds(M,t0=frame) 
                else:
                    vmin,vmax = check.bounds_pix(W)
        
            if Dt is not None:
                data = data/Dt
    
            if field in ['Ux','Uy']:
                vmax = np.abs(vmax)
                vmin = -np.abs(vmax)#*100
        
            if field in ['omega']:
                data = data
                vmax = np.abs(vmax)/5.#*15#/5.
                vmin = -np.abs(vmax)#*10#*100#vmax
            
            if field in ['strain']:
                data = data
                vmax = np.abs(vmax)/20.#*15#/5.
                vmin = -np.abs(vmax)#*10#*100#vmax
        
            if field in ['E']:
        #std = np.std(data[...,frame])
                vmax = vmax**2
                vmin = vmin**2
        
            if field in ['enstrophy']:    
                vmax = (vmax/5.)**2
                vmin = (vmin)**2
        
    if log:
        vmax = np.log10(vmax)
        if vmin>0:
            vmin = np.log10(vmin)
        else:
            vmin = vmax /100.
    n = (X.shape[0]-dimensions[0])/2
    if n!=0:
        X = X[n:-n,n:-n]
        Y = Y[n:-n,n:-n]        
    color_plot(X,Y,data[...,0],show=show,fignum=fignum,vmin=vmin,vmax=vmax,log=log,cbar=cbar) 
#    time_stamp(M,frame)
    if colorbar==True:
        plt.colorbar()
    
#    plt.axis('equal')
    if tstamp:
        t = M.t[frame]
        Dt = M.t[frame+1]-M.t[frame]
        s = ', t = '+str(np.round(t*1000)/1000)+' s, Dt = '+str(np.round(Dt*10000)/10)+'ms'
    else:
        s=''
        
    figs = {}
    figs.update(legende('X (mm)','Y (mm)',field+s,display=display,cplot=True,show=show))
    
    return figs
#    title = os.path.basename(M.Sdata.fileCine)

def movie(M,field,indices=None,compute=False,label='',Dirname = './',tracking=False,**kwargs):
    figs = {}
    if indices == None:
        nx,ny,nt = M.shape()
        indices = range(1,nt-1)
    if tracking:
        import stephane.vortex.track as track
        
    Dirname = Dirname + 'Movie_'+field+M.Id.get_id()+'/'
    
    fignum=1
    fig,ax = set_fig(fignum,subplot=111)
    plt.clf()
    
    start=True
    for frame in indices:
        figs.update(Mplot(M,field,frame,compute=compute,**kwargs))
        
        if tracking:
            tup = track.positions(M,frame,field='omega',display=False,sigma=3.,fignum=fignum)
        # print(tup)
#            graph([tup[0]],[tup[2]],label='bo',linewidth=3,fignum=fignum) #minimum position
            graph([tup[1]],[tup[3]],label='ro',linewidth=3,fignum=fignum) # maximum position
        
        
        if start:
            colorbar(label=label)
            start = False
        
        #print(Dirname)
        save_figs(figs,savedir=Dirname,suffix='_'+str(frame),dpi=100,frmt='png',display=False)
        
        plt.cla()
    
def time_stamp(M,i,x0=-80,y0=50,fignum=None):
    t = M.t[i]
    Dt = M.t[i+1]-M.t[i]
    s = 't = '+str(np.round(t*1000)/1000)+' s, Dt = '+str(np.round(Dt*10000)/10)+'ms'
 #   print(s)
    if fignum is not None:
        set_fig(fignum)
    plt.text(x0,y0,s,fontsize=20)
    
def vfield_plot(M,frame,fignum=1):
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
    x=M.x
    y=M.y
    Ux=M.Ux[:,:,frame]
    Uy=M.Uy[:,:,frame]
    
    colorCodeVectors=True
    refVector = 1.
    vectorScale = 100
    vectorColormap='jet'
    
    #bounds
    #chose bounds from the histograme of E values ?
    scalarMinValue=0
    scalarMaxValue=100
    
    # make the right figure active
    set_fig(fignum)
    
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
    
    legende('$x$ (mm)','$y$ (mm)','')
    
        # add reference vector
    #if settings.showReferenceVector:
    #        plt.quiverkey(Q, 0.05, 1.05, refVector, str(refVector) + ' m/s', color=settings.vectorColor)

    #overwrite existing colorplot
    refresh(False)

######################################################################
#################### Histograms and pdfs #############################
######################################################################

def hist(Y,Nvec=1,fignum=1,num=100,step=None,label='o-',log=False,normalize=True,xfactor=1,**kwargs):
    """
    Plot histogramm of Y values
    """
    set_fig(fignum)
   # print('Number of elements :'+str(len(Y)))
    if step is None:
        n,bins=np.histogram(np.asarray(Y),bins=num,**kwargs)
      #  print(bins)
    else:
        d = len(np.shape(Y))
#        print('Dimension : '+str(d))
        N=np.prod(np.shape(Y))
        if N<step:
            step=N
        n,bins=np.histogram(np.asarray(Y),bins=int(N/step))
        
    if normalize:
        dx = np.mean(np.diff(bins))
        n=n/(np.sum(n)*dx)
        
    xbin=(bins[0:-1]+bins[1:])/2/xfactor
    n = n*xfactor
       
    if log:
    # Plot in semilogy plot
        semilogy(xbin/Nvec,n,fignum,label)
    else:
        plt.plot(xbin,n,label)
        plt.axis([np.min(xbin),np.max(xbin),0,np.max(n)*1.1])
        
    refresh()
    return xbin,n    
    
def pdf(M,field,frame,Dt=10,Dx=1024,label='ko-',fignum=1,a=15.,norm=True,sign=1):
    import stephane.manager.access as access
    Up = access.get(M,field,frame,Dt=Dt)
    
    limits = [(0,Dx),(0,Dx)]
    Up = sign*access.get_cut(M,field,limits,frame,Dt=Dt)
    
    figs = distribution(Up,normfactor=1,a=a,label=label,fignum=fignum,norm=norm)
    
    return figs

def pdf_ensemble(Mlist,field,frame,Dt=10,Dx=1024,label='r-',fignum=1,a=10.,norm=True,model=False):

    import stephane.manager.access as access

    U_tot = []
    
    for M in Mlist:
        pdf(M,field,frame,Dt=Dt,Dx=Dx,label='k',fignum=fignum,a=a,norm=False)
        
        Up = access.get(M,field,frame,Dt=Dt)
       # limits = [(0,Dx),(0,Dx)]
    #    Up = access.get_cut(M,field,limits,frame,Dt=Dt) 
        # if Dx is larger than the box size, just keep all the data
        U_tot = U_tot + np.ndarray.tolist(Up)
        
    N = len(Mlist)
    U_tot = np.asarray(U_tot)
    
    x,y,figs = distribution(U_tot,normfactor=N,a=a,label=label,fignum=fignum,norm=norm)
    
    if model:
        n = len(y)
        b = y[n//2]    
        Dy = np.log((y[n//2+n//8] + y[n//2-n//8])/2./b)
    
        a = - Dy/x[n//2+n//8]**2
    
        P = b*np.exp(-a*x**2)
        semilogy(x,P,label='b.-',fignum=fignum)
    
    set_axis(min(x),max(x),1,max(y)*2)
    if field=='omega' or field=='strain':
        unit = ' (s^-1)'
    elif field=='E':
        unit = 'mm^2/s^2'
    else:
        unit = ' (mm/s)'
    figs = {}
    figs.update(legende(field+unit,field+' PDF',time_label(M,frame)))
    return figs
    
def avg_from_dict(dd,keyx,keyy,times,fignum=1,display=True,label='b-'):
    """
    Compute the average function from a dictionnary with keys (time,keyx) (time,keyy) for time in times
    
    """    
    avg = {}
    avg[keyx]= np.mean([dd[(time,keyx)] for time in times],axis=0)
    avg[keyy]= np.mean([dd[(time,keyy)] for time in times],axis=0)
    
    std = {}
    std[keyx]= np.std([dd[(time,keyx)] for time in times],axis=0)
    std[keyy]= np.std([dd[(time,keyy)] for time in times],axis=0)
    
    if display:
        for time in times:
            graph(dd[(time,keyx)],dd[(time,keyy)],label='k-',fignum=fignum,color='0.7')
                        
        errorbar(avg[keyx],avg[keyy],std[keyx],std[keyy],fignum=fignum,label=label)
    
    return avg,std

def distribution(Y,normfactor=1,a=10.,label='k',fignum=1,norm=True):

    Y = np.asarray(Y)
    Median = np.sqrt(np.nanmedian(Y**2))
    
    "test if the field is positive definite"
    t = Y>=0
    
    if norm:
        bound = a
    else:
        bound = a*Median
    step = bound/10**2.5

    if t.all():
        x = np.arange(0,bound,step)
    else:
        x = np.arange(-bound,bound,step)
    
    if norm:
        Y = Y/Median
        
    n,bins = np.histogram(Y,bins=x)
    xbin = (bins[:-1]+bins[1:])/2
    n = n/normfactor # in case of several series (ensemble average)

    semilogy(xbin,n,label=label,fignum=fignum)
    set_axis(min(xbin),max(xbin),min(n)/2,max(n)*2)
    figs = {}
    figs.update(legende('','PDF',''))
    
    val = 0.5
    x_center = xbin[np.abs(xbin)<val]
    n_center = n[np.abs(xbin)<val]
    
    moy = np.sum(n*xbin)/np.sum(n)
    std = np.sum(n*(xbin-moy)**2)/np.sum(n)
    
    print("Distribution : "+str(moy)+' +/- '+str(std))
    #    a = fitting.fit(fitting.parabola,x_center,n_center)
    #    n_th = fitting.parabola(xbin,a)
    #    graph(xbin,n_th,label='r-',fignum=fignum)
    return xbin,n,figs
####### to add :

#   1. streamline plots
#   2. Strain maps
#   3. Vorticity maps

#for i in range(10,5000,1):
#    vfield_plot(M_log[4],i,1)
   # input()
