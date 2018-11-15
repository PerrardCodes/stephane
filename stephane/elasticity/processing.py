



import glob
import PIL.Image as Image
import numpy as np

import scipy.ndimage.measurements as meas

import stephane.display.graphes as graphes
import stephane.display.panel as panel

def load(folder,step=1,substack=None):
    
    l = glob.glob(folder+'*')
    n = len(l)
    
    data = []
    if substack is None:
        substack=range(0,n,step)
    for i in substack:
        im = Image.open(l[i])
        pix = im.load()
        H,L = im.size
        data.append(np.asarray([[pix[i,j] for j in range(L)] for i in range(H)]))
        #print(im.size)
        
    D = {}
    D['im']=np.asarray(data)
    D['H']=H
    D['L']=L
    D['folder']=folder
    
    return D
    
def load_cine(c,i0=2,n=1100,step=100):
    data = []
    for i in range(i0,n,step):
        print(i)
        im = np.asarray(c.get_frame(i))
        H,L = im.shape
        data.append(im)
        #print(im.size)
        
    D = {}
    D['im']=np.asarray(data)
    D['H']=H
    D['L']=L
#    D['folder']=folder   
    return D
            
def interface(D,k,T=10):
    Hmin = 300
    Hmax = 1000
    
#    data = D['im'][k,slice(Hmin,Hmax),slice(D['L'])]
    data = D['im'][k,slice(D['L']),slice(Hmin,Hmax)]
    
    #print(data.shape)
    
    x = range(D['L'])
    y = range(D['H'])
    X,Y = np.meshgrid(x,y)
    
    grad = np.gradient(data)
    
    drop = np.sqrt(grad[0]**2+grad[1]**2)  #compute the norm of the gradient

    databin = drop>T #binarize the image

    label,num = meas.label(databin)
    sizes = [np.sum(label==i) for i in range(1,num)]
    indices = np.where(label==np.argmax(sizes)+1)
    return X[indices],Y[indices]

def profile(folder,i0=None,j0=400,step=1,end=None):
    
    l = glob.glob(folder+'*')
    n = len(l)
    
    if end is None:
        ilist = range(0,n,step)
    else:
        ilist = range(0,end,step)
        
    data = []
    for i in ilist:
        if i%1000==0:
            print(i)
        im = Image.open(l[i])
        pix = im.load()
        H,L = im.size
        if i0 is None:
            data.append(np.asarray([pix[i,j0] for i in range(H)]))
        else:
            data.append(np.asarray([pix[i0,j] for j in range(L)]))
        #print(im.size)
        
    D = {}
    D['im']=np.asarray(data)
    D['H']=H
    D['L']=L
    D['folder']=folder
    
    return D
    
def show(D,k):
    
    fig,axes = panel.make([111],fignum=1)
    fig.set_size_inches(D['L']/100,D['H']/100)

    x = range(D['L'])
    y = range(D['H'])
    X,Y = np.meshgrid(x,y)
    
    graphes.color_plot(X,Y,D['im'][k,...])
    #graphes.graph([400],[1200],label='ro')
    
#Xmean = np.mean(X*drop)/np.mean(drop)
#Ymean = np.mean(Y*drop)/np.mean(drop)
#print(Xmean,Ymean)
#N = 10**4
#T = [len(np.where(drop>k)[0]) for k in range(50,N)]
#graphes.graph(range(50,N),T)
#graphes.set_axis(0,N,0,20000)
#T = 100
    #indices = np.where(drop>T)
#indices_2 = np.where(label==np.argmax(sizes)+1)
    #graphes.graph(X[indices],Y[indices],label='k.',fignum=2)
#Xmean = np.mean(X[indices])#*drop[indices])/np.mean(drop[indices])
#Ymean = np.mean(Y[indices])#*drop[indices])/np.mean(drop[indices])
#graphes.graph([Xmean],[Ymean],label='ro',fignum=2)
