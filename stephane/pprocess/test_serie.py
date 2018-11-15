

import stephane.tools.rw_data as rw_data
import glob
import stephane.tools.browse as browse
import stephane.pprocess.check_piv as check
import stephane.mdata.Mdata_PIVlab as pivlab
import stephane.analysis.cdata as cdata
import stephane.display.graphes as graphes

import numpy as np
import stephane.analysis.vgradient as vgradient
import stephane.manager.access as access

def test_Dt(W):
    
    figs = {}
    ratio = []
    Dt_list = range(1,11)
    for i,Dt in enumerate(Dt_list):
        print('Dt : '+str(Dt))
        Dir = '/Users/stephane/Documents/Experiences_local/PIV_tests/Database/Turbulence/PIV_data/'
        base = '_2016_03_01_PIV_sv_vp_zoom_zoom_X10mm_M24mm_fps10000_n14000_beta500m_H1180mm_S150mm_1'
        Dirbase=Dir + 'PIVlab_ratio2_W'+str(W)+'pix_Dt_'+str(Dt)+base
        fileList = glob.glob(Dirbase+'/*.txt')        
        dataList = get_data(fileList)

        r,fig = test_bound(dataList,W,Dt,fignum=1+i)
        figs.update(fig)
        ratio.append(r)
    
    graphes.graph(Dt_list,ratio,fignum=11,label='ko')
    graphes.graph([0,11],[100,100],label='r-',fignum=11)
    graphes.set_axis(0,11,0,110.)
    figs.update(graphes.legende('Dt (# frame)','Percentage of good measurements',''))
        
def do(Mlist):
    """
    Run a serie of PIV tests on a given set of data
    """
    ratio = {}
    figs = {}
    for i,key in enumerate(Mlist.keys()):
        (Dt,W) = key 
        M = Mlist[key]
        r,fig = test_M(M,W,Dt,fignum=1+i)
        figs.update(fig)
        print(str((Dt,W))+': '+str(r))
        ratio[key] = r
        
    return ratio,figs
    
def accuracy(M,frames):
    """
    from a serie of adjacent frames compute :
        the ratio of measurements in the boundaries (Umin, Umax). 
    
        Compute the noise level on the velocity field by averaging over adjacent frames in time (hypothesis of well resolved dynamics)
        
        the ratio of measurements within the shear limit (dUmax)
        
        Compute the velocity gradient noise level using the same time-averaging technic
    """
    for frame in frames:
        Ux = access.get(M,'Ux',frame)
        
    
    
def get_data(fileList):    
    dataList = []
    for name in fileList:
        W = browse.get_number(name,'_W','pix_')
        Dt = int(browse.get_string(name,'pix_Dt_','_'))
        Header,Data=rw_data.read_dataFile(name,',',',')
    #    indexA=browse.get_number(Header[1],'A: im','.tiff') #name of the image is localized on the second line of the ASCII file
    #    indexB=browse.get_number(Header[1],'B: im','.tiff')
        Data = pivlab.switch_keys(Data)
        dataList.append(Data)
        for key in ['u','v']:
            Data[key] = (cdata.rm_nans([np.asarray(Data[key])],d=1,rate=0.05))[0]
    
    return dataList
    #print(Data.keys())
    #print(type(Data['u']))
    
    
def test_bound(dataList,W,Dt,**kwargs):
    maxn = 0
    Umin,Umax = bounds_pix(W)
    
    ratio = []
    for data in dataList:
#        values = np.asarray(data['u'])**2+np.asarray(data['v']**2)
        values = np.sqrt(np.asarray(data['u'])**2+np.asarray(data['v'])**2)
        r = len(np.where(np.logical_and(values>Umin,values<Umax))[0])*100./len(data['u'])
        ratio.append(r)
        xbin,n = graphes.hist(values,normalize=False,num=200,range=(0.,2*Umax),**kwargs)#xfactor = Dt
        maxn = max([maxn,max(n)*1.2])
        
    ratio = np.nanmean(np.asarray(ratio))
    graphes.graph([Umin,Umin],[0,maxn],label='r-',**kwargs)
    graphes.graph([Umax,Umax],[0,maxn],label='r-',**kwargs)
    graphes.set_axis(0,Umax*1.2,0,maxn)
    title = 'Dt = '+str(Dt)+', W = '+str(W)+'pix'
    fig = graphes.legende('U (pix)','Histogram of U',title)
   # graphes.set_axis(0,1.5,0,maxn)
    
    return ratio,fig
        #graphes.graph()

def test_M(M,W,Dt,**kwargs):
    maxn = 0
    Umin,Umax = check.bounds_pix(W)
    
    ratio = []
#        values = np.asarray(data['u'])**2+np.asarray(data['v']**2)
    values = np.sqrt(np.asarray(M.Ux)**2+np.asarray(M.Uy)**2)
    
    N = np.prod(values.shape)
    print(N)
    r = len(np.where(np.logical_and(values>Umin,values<Umax))[0])*100./N
    ratio.append(r)
    xbin,n = graphes.hist(values,normalize=False,num=200,range=(0.,2*Umax),**kwargs)#xfactor = Dt
    maxn = max([maxn,max(n)*1.2])
        
    ratio = np.nanmean(np.asarray(ratio))
    graphes.graph([Umin,Umin],[0,maxn],label='r-',**kwargs)
    graphes.graph([Umax,Umax],[0,maxn],label='r-',**kwargs)
    graphes.set_axis(0,Umax*1.2,0,maxn)
    title = 'Dt = '+str(Dt)+', W = '+str(W)+'pix'
    fig = graphes.legende('U (pix)','Histogram of U',title)
   # graphes.set_axis(0,1.5,0,maxn)
    
    return ratio,fig
    
        
def shear_limit_M(M,W,Dt,type=1,**kwargs):
    """
    Test the shear criterion : dU/W < 0.1 
    """
    values = access.get(M,'strain',frame)
    
    M,field = vgradient.compute(M,'strain',step=1,filter=False,Dt=1,rescale=False,type=type,compute=False)
    values = getattr(M,field)#/W
    
    dUmin,dUmax = check.shear_limit_M(M,W)
    
    xbin,n = graphes.hist(values,normalize=False,num=200,range=(-0.5,0.5),**kwargs)#xfactor = Dt
    maxn = max(n)*1.2
    
    graphes.graph([dUmin,dUmin],[0,maxn],label='r-',**kwargs)
    graphes.graph([dUmax,dUmax],[0,maxn],label='r-',**kwargs)
    graphes.legende('','','')
    
def example():
    pass
    
def main(W=16):
    test_Dt(W)
