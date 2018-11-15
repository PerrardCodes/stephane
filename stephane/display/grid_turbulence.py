

import stephane.display.graphes as graphes
import stephane.display.vfield as vfield
import stephane.analysis.vgradient as vgradient
import stephane.analysis.corr as corr

import stephane.analysis.length_scales as scale

import math
import numpy as np

import matplotlib.patches as patches
import matplotlib.pyplot as plt

import stephane.tools.browse as browse
import os.path

import stephane.analysis.Fourier as Fourier
import stephane.manager.access as access
#Codes for plotting data of grid turbulence experiments.
#include plots of box, grid, multiplot

def panel(M,fignum=1):
    fig=plt.figure(fignum)
    fig.set_size_inches(18.5, 16)
    
    #example of three subplots
    ax1 = box_axes(M,fig,131)   #for the side view
    graphes.legende('','','Side view')
    draw_fieldofview(M.Sdata,ax1,view='side')
    
    ax2 = box_axes(M,fig,132)   # for the front view
    graphes.legende('','','Front view')
    draw_fieldofview(M.Sdata,ax2,view='front')

    ax3 = box_axes(M,fig,133)   # for the front view
    graphes.legende('','','Front view')
    draw_fieldofview(M.Sdata,ax3,view='front')

#    ax3 = fig.add_subplot(233)
#    ax4 = fig.add_subplot(236)
    axes = [ax1,ax2,ax3]#,ax3,ax4]
    return axes
    
def panel_graphs(M,subplot=[2,2],fignum=1,W=24.5,L=16):
    fig=plt.figure(fignum)
    fig.set_size_inches(W, L)
    
    axes = []
    n = subplot[1]
    for j in range(subplot[1]):       
        for i in range(subplot[0]):
            num = subplot[0]*100+subplot[1]*10+i*n+j+1 
            axes.append(plot_axes(M,fig,num))
    
    return axes    
   
def flex_panel(M,fignum=1):
    fig=plt.figure(fignum)
    fig.set_size_inches(24.5, 16)
    
    #example of three subplots
    axes = []
    #ax1 = box_axes(M,fig,131)   #for the side view
    #graphes.legende('','','Side view')
    #draw_fieldofview(M.Sdata,ax1,view='side')
    #axes.append(ax1)
    #ax.append(plot_axes(M,fig,234))#3 = fig.add_subplot(236)
    
    ax2 = box_axes(M,fig,121)   # for the front view
    graphes.legende('','','Front view')
    if hasattr(M,'Sdata'):
        draw_fieldofview(M.Sdata,ax2,view='front')
    
    axes.append(ax2)
    
#    ax3 = box_axes(M,fig,143)   # for the front view
#    graphes.legende('','','Front view')
#    if hasattr(M,'Sdata'):
#        draw_fieldofview(M.Sdata,ax2,view='front')
#    axes.append(ax3)
    
    #Scales panel
    axes.append(plot_axes(M,fig,222))
    axes.append(plot_axes(M,fig,224))#3 = fig.add_subplot(236)
    
    #Distribution panel
#    axes.append(plot_axes(M,fig,233))#3 = fig.add_subplot(236)    
#    axes.append(plot_axes(M,fig,236))
    
    #3 = fig.add_subplot(236)
#    axes = [ax1,ax2,ax3,ax4]#,ax3,ax4]
    return axes 
    

def plot_axes(M,fig,num):
    ax = fig.add_subplot(num)
  #  ax.set_aspect('equal', adjustable='box')
    #graphes.legende('','','Front view')
    #draw_fieldofview(M.Sdata,ax3,view='front')
    return ax

def box_axes(M,fig,num):
    
    ax = plot_axes(M,fig,num)
    ax.axis('off')
    draw_box(M,ax)
    
    if hasattr(M,'param'):
        S = -M.param.stroke
        draw_grid(M,ax,0,0,facecolor='w')
        draw_grid(M,ax,0,S)
    
    return ax
    
def draw_box(M,ax):
    """
    Draw the box. Dimensions are given by the date
    """
    offset = 150.
    width = 250.
    height = 700.
    
    eps = width/100
    ax.set_xlim([-width/2-eps,width/2+eps])
    ax.set_ylim([offset-height-eps,offset+eps])
    ax.set_ylim([offset-height-eps,50])
   
    c = tuple([0.8 for i in range(3)])     
   
    ax.add_patch(patches.Rectangle((-width/2,-height+offset),width,height,fill=False,linewidth=10,edgecolor=c))
    #plt.draw()
    
def draw_fieldofview(S,ax,fignum=None,color='r',view='side'):
    """
    Draw the field view of the camera, currently side view only.
    """
    if fignum is not None:
        graphes.set_fig(fignum)
  #  print(S.fileCine)

    u0,L=field_of_view(S)
    if view=='side':
        graphic_run(ax,u0[1],u0[2],L[1],L[2],color=color)
    if view=='front':
        graphic_run(ax,u0[0],u0[2],L[1],L[2],color=color)
        graphic_run(ax,u0[0]+L[0],u0[2],L[1],L[2],color=color)
        graphic_run(ax,u0[0],u0[2],L[0],L[1],color=color)
        graphic_run(ax,u0[0],u0[2]+L[2],L[0],L[1],color=color)
    if view=='bottom':
        pass
        
        
    #find the plane of measurement from Sdata file and     
        
def field_of_view(S):
    """
    From a Sdata, return the field of view given by : 
    down left corner
    length,width,height
    """
    if hasattr(S.param,'im_ref'):
        angle = S.param.angle*math.pi/180
        fx = S.param.fx
        
        (Ypix,Xpix) = np.shape(S.param.im_ref)
        
        #find length,width and depth
        if S.param.typeview=='sv':
            X = (Ypix*math.sin(angle)+Xpix*math.cos(angle))*fx #along x
            Y = 0    #along y
            Z = (Xpix*math.sin(angle)+Ypix*math.cos(angle))*fx #along Z

        if S.param.typeview=='bv':
            #no rotation needed (X and Y are considered as equivalent directions)
            X = (Ypix*math.sin(angle)+Xpix*math.cos(angle))*fx
            Y = (Xpix*math.sin(angle)+Ypix*math.cos(angle))*fx
            Z = 0
            
        #find 0 position (lower left corner)
        if S.param.typeview == 'sv':
            #vertical plane, the horizontal position is given py Xplane
            x0 = (S.param.x0*math.cos(angle)+S.param.y0*math.sin(angle))*fx-X
            y0 = S.param.Xplane
            z0 = -(S.param.x0*math.sin(angle)+S.param.y0*math.cos(angle))*fx

        if S.param.typeview == 'bv':
            #vertical plane, the horizontal position is given py Xplane
            x0 = (S.param.x0*math.cos(angle)+S.param.y0*math.sin(angle))*fx-X
            y0 = -(S.param.x0*math.sin(angle)+S.param.y0*math.cos(angle))*fx
            z0 = S.param.Zplane
            
        u0 = tuple([x0,y0,z0])
        L = tuple([X,Y,Z])
        
        return u0,L
    else:
        print('No reference image given')
        return None        

def graphic_run(ax,x0=0,z0=0,X=0,Z=200,color='r'):
    #a run is represented graphycally by a red line with two small black line at each extremity

#    style = patch.ArrowStyle.BarAB()
    arrow = patches.Arrow(x0,z0,X,Z,width=.0,color=color,linewidth=1)
    ax.add_patch(arrow)
    
    ratio=20        
    ledge = patches.Arrow(x0-Z/ratio,z0-X/ratio,2*Z/ratio,2*X/ratio,width=.0,color='k',linewidth=1)
    redge = patches.Arrow(x0+X-Z/ratio,z0+Z-X/ratio,2*Z/ratio,2*X/ratio,width=.0,color='k',linewidth=1)
    ax.add_patch(ledge)
    ax.add_patch(redge)
    
    graphes.refresh()    
    

def draw_grid(M,ax,x0,y0,Mesh=50,Bar=10,H=10,N=5,facecolor='y',edgecolor='k',linewidth=2):
    """
    Draw the grid at the specified position.
    """
    file = M.Sdata.fileCine
    ind = browse.get_string(os.path.basename(file),'_M',end='mm')
    
    if M.Id.date=='2016_03_05':
        ind = 24
        
    if not ind=='':
        Bar = 4
        H = 8
        Mesh = int(ind)#print(int(ind))
        N = 9
    else:
        Mesh = 50
        
    positions = [(i*Mesh+x0,y0) for i in range(-(N//2),N//2+1)]
    width = Bar
    height = H
   # grid = []
    for p in positions:
        square = patches.Rectangle(centered(p,width,height),width,height,facecolor=facecolor,edgecolor=edgecolor,linewidth=linewidth)    
        ax.add_patch(square)
        
def centered(t,width,height):
    return (t[0]-width/2,t[1]-height/2)
    
def make_movie(M,Range=None,field=['E','vorticity'],Dirbase=None,Dt=1):
    if Dirbase==None:
        Dirbase = '/Users/stephane/Documents/Experiences_local/Accelerated_grid/PIV_data' #local saving
    
    nx,ny,nt = M.shape()
    
    axes = panel(M)
    if Range is not None:
        start,stop,step = tuple(Range)
    else:
        start,stop,step = tuple([0,nt,1])
    frames = range(start,stop,step)
    
#    field = ['Ux','Uy']
  #  field = ['E','vorticity']
    figs={}
    Dirname = Dirbase+'/'+M.Id.get_id()+'/'+graphes.remove_special_chars(str(field))+'/'
    print(Dirname)
    
    for i,f in enumerate(field):
        if not hasattr(M,f):
            M,field_n[i] = vgradient.compute(M,field[i],Dt=Dt)
            field = field_n
            
#    M,field[1] = vgradient.compute(M,field[1],Dt=Dt)
    
    if Dt>1:
        print('Smoothed data')
        Dirname = Dirname + 'Smooth_Dt_'+str(int(Dt))+'/'

    for i in frames:        
        graphes.set_fig(1)
        for axe in axes:
            plt.sca(axe)
            plt.cla()
            
        axes = panel(M,fignum=1)        
        plt.sca(axes[1])
        print(i)
        graphes.Mplot(M,field[0],i,fignum=1,log=True)
        plt.text(-10,80,field[0],fontsize=20)

        plt.sca(axes[2])
        graphes.Mplot(M,field[1],i,fignum=1,log=True)
        plt.text(-10,80,field[1],fontsize=20)
        
        figs.update(graphes.legende('','','Front view',cplot=True))
        graphes.save_figs(figs,savedir=Dirname,suffix='_'+str(i),dpi=100,frmt='png',display=False)
    return axes
    
def subplot(ax,x,y_list,labels=None,fignum=1):    
    pass
    
def select_range(M,Range):
    nx,ny,nt = M.shape()
    
    if Range is not None:
        start,stop,step = tuple(Range)
    else:
        start,stop,step = tuple([0,nt,1])
    frames = range(start,stop,step)
    
    return frames
    
def make_plot(M,Range=None,color='k',field=['E','vorticity'],Dirbase=None,Dt=1,example=False,total=True,fignum=1,save=True):
    if Dirbase==None:
        Dirbase = '/Users/stephane/Documents/Experiences_local/Accelerated_grid/PIV_data/Test6/' #local saving
        Dirbase = './Stat_avg/Panel/'+M.Id.date
    
    axes = flex_panel(M,fignum=fignum)
    
   # for axe in axes:
    #    plt.sca(axe)
        #plt.cla()
    
    frames = select_range(M,Range)
    
#    field = ['Ux','Uy']
  #  field = ['E','vorticity']
    figs={}
    
    if hasattr(M,'id'):
        Dirname = Dirbase+'/'+M.Id.get_id()+'/'+graphes.remove_special_chars(str(field))+'/'
    else:
        Dirname = Dirbase+'/JHTD_Data/'+graphes.remove_special_chars(str(field))+'/'
    print(Dirname)
    
    if Dt>1:
        print('Smoothed data')
        Dirname = Dirname + 'Smooth_Dt_'+str(int(Dt))+'/'

    #Dt = 50
    t_moy,E_moy = access.avg_vs_t(M,'E',frames,Dt=Dt)
    t_moy,Omega_moy = access.avg_vs_t(M,'omega',frames,Dt=Dt)
    t_moy,Strain_moy = access.avg_vs_t(M,'strain',frames,Dt=Dt)
    t_moy,Y_pdf_moy = access.avg_vs_t(M,field[1],frames,Dt=Dt)
    
    epsilon = scale.dissipation_rate(Omega_moy) #dissipation rate
    eta = scale.K_scale(epsilon)  #dissipative scale
    
    micro_1 = np.sqrt(1./4*E_moy/Strain_moy**2)
    micro_2 = np.sqrt(E_moy/(Omega_moy)**2) 
    
    Re = scale.Re(micro_1,eta)
    Re_lambda_1 = scale.Re_lambda(E_moy,micro_1)
    Re_lambda_2 = scale.Re_lambda(E_moy,micro_2)
    
    L = scale.I_scale(Re,E_moy)
    
#    plt.sca(axes[2])
#    graphes.graphloglog(t_moy,E_moy,fignum=1,label=color+'o-')
    #Fourier.add_theory(t_moy,E_moy,[-2.],fignum=1)
#    graphes.graphloglog(t_moy,Y_pdf_moy,fignum=1,label=color+'s-')
#    graphes.graphloglog(t_moy,epsilon,fignum=1,label='gv-')
    #Fourier.add_theory(t_moy,epsilon,[-3.],fignum=1)    
#    figs.update(graphes.legende('t (s)','E (mm^2/s^2), epsilon(mm^2/s^-3)',''))
    
    plt.sca(axes[1])
    graphes.graphloglog(t_moy,eta,fignum=fignum,label=color+'o-')
    graphes.graphloglog(t_moy,micro_1,fignum=fignum,label=color+'s-')
    graphes.graphloglog(t_moy,micro_2,fignum=fignum,label='cp-')
    graphes.graphloglog(t_moy,L,fignum=fignum,label='gv-')
    figs.update(graphes.legende('t (s)','eta (mm), lambda (mm)',''))        
    
    plt.sca(axes[2])
    graphes.graphloglog(t_moy,Re,fignum=fignum,label=color+'o-')
    graphes.graphloglog(t_moy,Re_lambda_1,fignum=fignum,label=color+'s-')
    graphes.graphloglog(t_moy,Re_lambda_2,fignum=fignum,label='cp-')
    figs.update(graphes.legende('t (s)','Re , Re_lambda',''))
    
   # print(t_moy)
#    print(Y_moy)
   # print(t)
#    print(Y_moy)
    indices = [0,1]
#    indices = [1,4]
    cla_axes = [axes[i] for i in indices]
    
    if save:
        graphes.save_figs(figs,savedir=Dirname,prefix='General',suffix='_vs_t',dpi=300,frmt='png',display=True)

    individual=False
    
    if example:
        frames_disp =[1200]
    else:
        step =frames[1]-frames[0]
        frames_disp = range(frames[0]+step*10,frames[-1],step*10)
        
    if individual:
        for frame in frames_disp:  
            #print(frame)
            #print(frames)
            i = frames.index(frame)      
            graphes.set_fig(1)
            for axe in cla_axes:
                plt.sca(axe)
                plt.cla()
            axes = flex_panel(M,fignum=1)   
            
            
            if total:     
                plt.sca(axes[0])
                graphes.Mplot(M,field[0],frame,fignum=1,log=True)
                plt.text(-10,80,field[0],fontsize=20)
                
            plt.sca(axes[1])
         #   graphes.graphloglog(t_moy,eta,fignum=1,label='ko-')
         #   graphes.graphloglog(t_moy,micro,fignum=1,label='ks-')
            graphes.graph([t_moy[i]],[eta[i]],fignum=1,label='ro')
            graphes.graph([t_moy[i]],[micro_1[i]],fignum=1,label='rs')
            graphes.graph([t_moy[i]],[micro_2[i]],fignum=1,label='rp')
            graphes.graphloglog([t_moy[i]],[L[i]],fignum=1,label='rv-')
            
            plt.sca(axes[2])
            graphes.graphloglog([t_moy[i]],[Re[i]],fignum=1,label='ro-')
            graphes.graphloglog([t_moy[i]],[Re_lambda_1[i]],fignum=1,label='rs-')
            graphes.graphloglog([t_moy[i]],[Re_lambda_2[i]],fignum=1,label='rp-')

        
            plt.sca(axes[3])
            figs.update(graphes.pdf(M,'Ux',frame,Dt=Dt,fignum=1,label='m^'))
            figs.update(graphes.pdf(M,'Uy',frame,Dt=Dt,fignum=1,label='b>'))
            figs.update(graphes.legende('t (s)','Re , Re_lambda',''))

            plt.sca(axes[4])
            figs.update(graphes.pdf(M,field[1],frame,Dt=Dt,fignum=1,label=color+'-'))
        
#        graphes.Mplot(M,field[1],i,fignum=1,log=True)
#        plt.text(-10,80,'PDF '+field[1],fontsize=20)       
#        figs.update(graphes.legende('','','Front view',cplot=True))
            graphes.save_figs(figs,savedir=Dirname,suffix='_'+str(frame),dpi=100,frmt='png',display=False)
            
    return axes
    
def make_plot_lin(Mlist,Range=None,color='k',label=None,field=[['Ux','Uy'],['omega']],Dirbase=None,Dt=1,example=False,total=True,fignum=1,save=True):
    M = Mlist[0]
    if Dirbase==None:
        Dirbase = '/Users/stephane/Documents/Experiences_local/Accelerated_grid/PIV_data/Test6/' #local saving
        Dirbase = './Stat_avg/Panel/'+M.Id.date
    
    axes = panel_graphs(M,subplot=[2,3],fignum=fignum)

    frames = select_range(M,Range)
    
    figs={}    
    if hasattr(M,'Id'):
        Dirname = Dirbase+'/'+M.Id.get_id()+'/'+graphes.remove_special_chars(str(field))+'/'
    else:
        Dirname = Dirbase+'/JHTD_Data/'+graphes.remove_special_chars(str(field))+'/'
    print(Dirname)
    
    if Dt>1:
        print('Smoothed data')
        Dirname = Dirname + 'Smooth_Dt_'+str(int(Dt))+'/'

    figs.update(plot_scales(Mlist,axes,fignum=fignum,color=color,label=label))
    
    plt.sca(axes[2])
    frame = 1500
    Dt = 1400
    
    if label is None:
        labels = ['m^','b>','ko']
    else:
        labels = [label,label,label]
        
    for i,f in enumerate(field[0]):#should contain either one or two fields
        figs.update(graphes.pdf_ensemble(Mlist,f,frame,Dt=Dt,fignum=fignum,label=labels[i],norm=False))
        figs.update(graphes.legende(f,'pdf of '+f,''))

    plt.sca(axes[3])
    for f in field[1]:
        figs.update(graphes.pdf_ensemble(Mlist,f,frame,Dt=Dt,fignum=fignum,label=labels[2],norm=False))
        figs.update(graphes.legende(f,'pdf of '+f,''))    

    plt.sca(axes[4])
    corr.corr_v_t(Mlist,frame,axes=['Ux','Ux'],N=200,p=1,display=True,save=False,label=labels[0],fignum=fignum)
    corr.corr_v_t(Mlist,frame,axes=['Uy','Uy'],N=200,p=1,display=True,save=False,label=labels[1],fignum=fignum)
    
    plt.sca(axes[5])
    corr.corr_v_t(Mlist,frame,axes=['omega','omega'],N=200,p=1,display=True,save=False,label=labels[2],fignum=fignum)
    
    if save:
        graphes.save_figs(figs,savedir=Dirname,prefix='General',suffix='_vs_t',dpi=300,frmt='png',display=True)
    else:
        return figs,Dirname

def make_stat_avg(Mlist,**kwargs):
    
    t_moy,eta,micro,L,Re_lambda,Re,epsilon = scales(Mlist)
    
    
def plot_scales(Mlist,axes,fignum=1,color='k',label=None):
    figs = {}
    
    t_moy,eta,micro,L,Re_lambda,Re,epsilon = scales(Mlist)
    
    if label is None:
        labels = [color+'o-',color+'s-','gv-']
    else:
        labels = [label,label,label]
    
    plt.sca(axes[0])
    graphes.semilogy(t_moy,eta,fignum=fignum,label=labels[0])
    graphes.semilogy(t_moy,micro,fignum=fignum,label=labels[1])
    # graphes.graphloglog(t_moy,micro_2,fignum=fignum,label='cp-')
    graphes.semilogy(t_moy,L,fignum=fignum,label=labels[2])
    figs.update(graphes.legende('t (s)','eta (mm), lambda (mm)',''))        
    
    plt.sca(axes[1])
    graphes.semilogy(t_moy,Re,fignum=fignum,label=labels[0])
    graphes.semilogy(t_moy,Re_lambda,fignum=fignum,label=labels[1])
     #graphes.graphloglog(t_moy,Re_lambda_2,fignum=fignum,label='cp-')
    figs.update(graphes.legende('t (s)','Re , Re_lambda',''))
    
    return figs
     
     
def make_stat_panel(Mlist,fignum=1,**kwargs):
    
    for M in Mlist:
        make_plot_lin([M],fignum=fignum,**kwargs)
    
    figs,Dirname = make_plot_lin(Mlist,fignum=fignum,label='r',**kwargs)
    graphes.save_figs(figs,savedir=Dirname,prefix='General',suffix='_vs_t',dpi=300,frmt='png',display=True)
    
def scales(Mlist,Range=None,Dt=1):
    """
    From a measurement set, return the three relevant lengths as a function of time.
    """
    M = Mlist[0]
    frames = select_range(M,Range)
    
    n = len(Mlist)
    E = [None for i in range(n)]
    Omega = [None for i in range(n)]
    Strain = [None for i in range(n)]
    
    for i,M in enumerate(Mlist):
    #Dt = 50
        t_moy,E[i] = access.avg_vs_t(M,'E',frames,Dt=Dt)
        t_moy,Omega[i] = access.avg_vs_t(M,'omega',frames,Dt=Dt)
        t_moy,Strain[i] = access.avg_vs_t(M,'strain',frames,Dt=Dt)
    E_moy = np.nanmean(np.asarray(E),axis=0)
    Omega_moy = np.sqrt(np.nanmean(np.power(np.asarray(Omega),2),axis=0))
    Strain_moy = np.sqrt(np.nanmean(np.power(np.asarray(Strain),2),axis=0))
    
    epsilon = scale.dissipation_rate(Omega_moy) #dissipation rate
    
    eta = scale.K_scale(epsilon)  #dissipative scale
    micro_1 = np.sqrt(E_moy/Strain_moy**2) #Taylor microscale

    Re = scale.Re(micro_1,eta)  #Reynolds number
    Re_lambda_1 = scale.Re_lambda(E_moy,micro_1)  #Taylor Reynolds number
 
    L = scale.I_scale(Re,E_moy) # integral length scale
    
    return t_moy,eta,micro_1,L,Re_lambda_1,Re,epsilon
    

def scale_select(M,selection,Range=None,Dt=1):
    table = {'eta':1,'micro':2,'L':3,'Re_lambda':4,'Re':5,'dissipation':6}
    
    tup = scales(M,Range=Range,Dt=Dt)
    
    if selection in table.keys():
        return tup[0],tup[table[selection]]
    else:
        print("Scale not found")
        return None
    #micro_2 = np.sqrt(E_moy/(Omega_moy)**2) 
    

def example(M,Range=None,field=['E','vorticity']):
    nx,ny,nt = M.shape()
    axes = panel(M)
    if Range is not None:
        start,stop,step = tuple(Range)
    else:
        start,stop,step = tuple([0,nt,1])
    frames = range(start,stop,step)
    
#    field = ['Ux','Uy']
  #  field = ['E','vorticity']
    
    figs={}
    
    Dirbase = '/Users/stephane/Documents/Experiences_local/Accelerated_grid/PIV_data' #local saving
    Dirname = Dirbase+'/'+M.Id.get_id()+'/'+graphes.remove_special_chars(str(field))+'/'
    print(Dirname)
    
    M,field[0] = vgradient.compute(M,field[0])
    M,field[1] = vgradient.compute(M,field[1])
    
    for i in frames:        
        graphes.set_fig(1)
        for axe in axes:
            plt.sca(axe)
            plt.cla()
        axes = panel(M,fignum=1)    
            
        plt.sca(axes[1])
        graphes.Mplot(M,field[0],i,log=True)
        plt.text(-10,20,field[0],fontsize=20)

        plt.sca(axes[2])
        graphes.Mplot(M,field[1],i,fignum=1,log=True)
        plt.text(-10,20,field[1],fontsize=20)
        
        figs.update(graphes.legende('','','Front view',cplot=True))
        graphes.save_figs(figs,savedir=Dirname,suffix='_'+str(i),dpi=100,frmt='png')
    return axes
    
    