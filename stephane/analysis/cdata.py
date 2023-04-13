# -*- coding: utf-8 -*-
"""
Created on Thu Jul 16 11:59:31 2015

@author: stephane
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Jul 16 11:52:45 2015

@author: stephane
"""

import numpy as np
import scipy.interpolate
import math
#import stephane.display.graphes as graphes
import scipy.ndimage.filters as filters
import time

### some functions to smooth results, remove NaN values, aberrant values and so on,
### based on adjacent values (useful for 2d arrays or greater dimension)


def nancount(data):
    print(data.shape)
    Nnan = np.sum(np.ndarray.flatten(np.isnan(data)))
    N = np.prod(data.shape)
    print('Ratio nan :' +str(Nnan/N))
    print('nan number :'+str(Nnan))
### CHarge un fichier Data## .
### CHarge un fichier Data## .

def remove_nan_3d(ff):
    #Fields with d components along the last axis
    Nc = ff.shape[-1]
    Nt = ff.shape[0]
    for i in range(Nt):
        for j in range(Nc):
            data = np.squeeze(ff[i,...,j])
            indices = np.where(np.isnan(data))
            #adapted to 3d data
            (t0,t1,t2) = indices
            for tup in zip(t0,t1,t2):
                data = replace_nan(data,tup)
            ff[i,...,j]=data
    return ff
        
def rm_nan(M,field,rate=0.1):
    """
    Remove nan values from a given field of M (one of its attribute)
    INPUT
    -----
    M : Mdata object
        dataset
    field : string
        attribute of M containing the data to be processed
    OUTPUT
    -----
    M : Mdata object
        dataset with removes nan values
    """
    U = getattr(M,field)
    Ulist =rm_nans([U],rate=rate)
    setattr(M,field,Ulist[0])

    return M

def rm_nans(U,d=2,rate=0.01):
    """
    look for NaN values and replace by an average of adjacent values
    INPUT
    -----
    U : list of np array
    OUTPUT
    -----
    U_new : list of np array
        nan removed
    """
    #
    U_new=[]
    for u in U:
        tlist=np.where(np.isnan(u))
        dimension=len(tlist)

        if dimension==0 or len(tlist[0])==0:
            return U
        if dimension==1:
            for t in tlist[0]:
             #   print(t)
                u=replace_nan(u,[t])
                U_new.append(u)
        if dimension==d+1:
            (i,j,k)=tlist

            nx,ny,nt=u.shape
            ratio=len(i)*1./(nx*ny*nt)
          #  print('Percentage of values to be removed : '+str(ratio*100))
            if ratio<rate:
                for p,t in enumerate(zip(i,j,k)):
            #        if np.mod(p,len(i)/10)==0:
#                        print(str(p*100/len(i))+ ' %')
                    u=replace_nan(u,t)
            else:
                print('Too many values to correct. Skip')
            U_new.append(u)

    #give the percentage or removed values
    return U_new


def replace_nan(u,t,mode='mean',b=1):
    """
    replace a NaN value by an average of non nan adjacent values (weighted ??)
    INPUT
    -----
    u : np array
        data to be nan removed
    t : tuple
         tuple of index positions
    mode : string
         'mean' or 'median'
    OUTPUT
    -----
    u : nan removed np array
    """
    #
    # if on a corner , mirror the other side ??
    # u corresponds to the matrix given
    # t corresponds to the indices of nan value to remove in u
    neigh,ind=fast_neighboors(u,t,b=b)
    if mode=='mean':
        u[t]=np.nanmean(neigh)
    if mode=='median':
        u[t]=np.nanmedian(neigh)
    return u
    #ind = where(logical_not(np.isnan(y)))[0]
#y1 = interp(range(len(y)), ind, y[ind])
#y1 = y1[ind[0]:ind[-1]]
def neighboors(u,t,b=1,bmin=0,return_d=False,t_dilate=1):
    """

    """
    #case b=1. for large b, corner and edge might be increased
    #corner = (t in [(0,0),(ny-1,0),(0,nx-1),(ny-1,nx-1)])
    #edge = ((i in [0,nx-1]) or (j in [0,ny-1])) and (not corner)
    #bulk = not (corner or edge)
    #only d=2 case. might be extended to higher dimension (up to 4 !)
    d=2
    dimensions=u.shape
#    neigh=u[j-b:j+b+1,i-b:i+b+1]
        #number of neighbours per line
#    n=2*b+1
#    neigh_vec=np.reshape(neigh,(n**2,)) # 8 neighbours + the central point
#    neigh_vec0=np.concatenate((neigh_vec[:n*b+b],neigh_vec[n*b+b+1:])) # list of the 8 neighbouhrs elements
    if len(dimensions)>=2:
      #  print(t)
        boolean=[[(t[0]-q)**2+(t[1]-p)**2<=b*d and (t[0]-q)**2+(t[1]-p)**2>bmin  for p in range(dimensions[1])] for q in range(dimensions[0])]
        ind = np.where(boolean)
      #  nan_neigh=np.where(np.isnan(ind))
      #  print('Number of neighboors : '+str(len(ind[0])))
      #  print('Number of NaN neighboors : '+str(len(nan_neigh[0])))
        if len(dimensions)==3:
            ind=ind+(t[2]*np.ones((len(ind[0]),),dtype=int),)        #might be useful to implement high order dimensions (4d)
        #d can be used later to weight the average by closest point
        if return_d:
            d=[(t[0]-tup[0])**2+(t[1]-tup[1])**2 for tup in ind]
            return u[ind],ind,d
        else:
         #   print(ind)
        #    print(ind)
            return u[ind],ind

    if len(dimensions)==1:
        ind = [t+i for i in range(-b,b+1)]
        ind.remove(t)

        return u[ind],ind

def fast_neighboors(u,t,b=1,bmin=0,return_d=False,t_dilate=1):
    """
    Find the neighboors of a given point in 2 or 3 dimensions
    """
    #case b=1. for large b, corner and edge might be increased
    #corner = (t in [(0,0),(ny-1,0),(0,nx-1),(ny-1,nx-1)])
    #edge = ((i in [0,nx-1]) or (j in [0,ny-1])) and (not corner)
    #bulk = not (corner or edge)

    dim = u.shape
    d = len(dim)
   # print('Dimension : '+str(d))
#    (i0,j0) = (t[0],t[1])

    bmin=[]
    bmax=[]
    for p in range(d):
        bmin.append(max([0,t[p]-b]))
        bmax.append(min([dim[p],t[p]+b+1]))

    if d==1:
        tab = [i for i in range(bmin[0],bmax[0])]
        tab.remove(t[0])
    if d==2:
        tab = [(i,j) for i in range(bmin[0],bmax[0]) for j in range(bmin[1],bmax[1])]
        tab.remove(t)
    if d==3:
        tab = [(i,j,k) for i in range(bmin[0],bmax[0]) for j in range(bmin[1],bmax[1]) for k in range(bmin[2],bmax[2])]
        tab.remove(t)

    ind = np.asarray(tab)#np.reshape(np.asarray(tab),((xmax-xmin)*(ymax-ymin),2))
    U_sec = np.asarray([u[t] for t in tab])
#     ymin = max([0,j0-b])
#    ymax = min([ny,j0+b+1])
#    central = np.where(ind==[t[0],t[1]])
    return U_sec,ind

def rm_data(M,field):
    if field=='velocity':
        return rm_data_2d(M)
    else:
        U = getattr(M,field)
        nx,ny,nt = U.shape

        val = 5
        for t in range(nt):
            for i in range(nx):
                for j in range(ny):
                    u,ind = neighboors(U,(i,j),b=1.5)
                    moy_u =np.nanmean(u)
                    std_u = np.nanstd(u)
           # print(str(moy_u)+' +/- '+str(std_u))
                    if (U[i,j]-moy_u)>(val*std_u):
                        U[i,j] = moy_u
                        print((U[i,j]-moy_u)/moy_u)

def local_smooth(M,t,Dt,borne=2.):
    """
    Replace aberrant values by the mean of their neighboors.
    Aberrant values are detected with test_neigh function
    INPUT
    -----
    M : Mdata object
        Must contain attributes Ux, Uy and a methode shape()
    t : int
        time index
    Dt : int
        number of frames to process
    OUTPUT
    -----
    None. Directly modify the M object given in argument
    """
    M.get('E')
    tup,error = rm_baddata(M.E[...,t:t+Dt],d_max=2,b=1,borne=borne,display=False)

    #tup_list_y,error_y = rm_baddata(M.Uy[...,t:t+Dt],d_max=2,b=1,borne=borne,display=False)
    M.Ux[...,t:t+Dt] = replace_data(M.Ux[...,t:t+Dt],tup)
    M.Uy[...,t:t+Dt] = replace_data(M.Uy[...,t:t+Dt],tup)
    return M

def remove_spurious(M,borne=2.,display=False):
    M.get('E')
    nx,ny,nt = M.shape()
    i = nt/2

    #only a small data chunk
    #tup,error = rm_baddata(M.E[...,i-5:i+5],d_max=2,b=1,borne=borne,display=False)
    #for j,t in enumerate(tup):
    #    tup[j] = (t[0],t[1],t[2] + i-5)
    tup,error = rm_baddata(M.E,d_max=2,b=1,borne=borne,display=False)

    if display:
        import stephane.display.graphes as graphes
        graphes.Mplot(M,'E',i,colorbar=True)
        xc,yc = [],[]
        for t in tup:
            xc.append([M.x[t[0],t[1]]])
            yc.append([M.y[t[0],t[1]]])
        graphes.graph(xc,yc,label='ro')

    #tup_list_y,error_y = rm_baddata(M.Uy[...,t:t+Dt],d_max=2,b=1,borne=borne,display=False)
    M.Ux = replace_data(M.Ux,tup)
    M.Uy = replace_data(M.Uy,tup)
    return M

#def rm_baddata(M,field,t,Dt,d_max=2,b=1,borne=5):
    #A = getattr(M,field)
    #U = A[...,t:t+Dt]


def replace_data(data,tup_list):
    """
    Replace a data point by the mean of its neighboors values
    """
    dim = data.shape
    if len(dim)==2:
        data = np.transpose(np.tile(data,(1,1,1)),(1,2,0))

    for tup in tup_list:
       # print(tup)
        data_neigh = fast_neigh(data,tup,d_max=2,b=1,rm_center=True)
        data[tup] = np.nanmean(data_neigh)
    return data

def rm_baddata(U,d_max=2,b=1,borne=1.2,display=True):
    """
    Local test (based on the neighboorhood of each point)
    """
    t0 = time.time()
    dim = U.shape
    if len(dim)==2:
        U = np.transpose(np.tile(U,(1,1,1)),(1,2,0))
        dim = U.shape
    error1 = np.zeros(dim)

    c=0.
    tup_list = []
    for i in range(dim[0]):
        for j in range(dim[1]):
            for k in range(dim[2]):
                tup = (i,j,k)
                errors = test_neigh(U,tup,d_max=d_max,b=b)
                if np.abs(errors)>borne:
                    c+=1
                    tup_list = tup_list + [tup]
                    error1[i,j,k]=0.
                else:
                    error1[i,j,k]=errors
    t1 = time.time()

    print('Percentage of aberrant values : '+str(c/np.prod(dim)*100) + ' %')
    print('Time elapsed per frame : '+str((t1-t0)/dim[2]) +' s')

    #if display:
    #   graphes.hist(error1,num=10**3,fignum=1,log=True)
 #   graphes.set_axis(-5,5,0,5)
    #graphes.hist(error2,num=10**3,fignum=2)
    return tup_list,np.nanstd(error1)
    #print(np.nanmedian(error2))


def test_neigh(data,tup,d_max=None,b=1):
    """
    Compute the vector averaged of the neighbourhood of data[t], and compare to its value.
    If it is far from the median value, replace by the mean vector
    INPUT
    -----
    data : np array
        full data set ? or only slice ? -> only slice, ie containing dimensions on which we want to average
        last dimension then corresponds to vector (if it is a vector quantity)
    t : tuple
        index of the center point

    OUTPUT
    -----
    """
    bound = 5

    data_neigh = fast_neigh(data,tup,d_max=d_max,rm_center=True,b=b)
    val = data[tup]
    moy = np.nanmedian(data_neigh)
    std = np.nanstd(data_neigh)

    error = (moy-val)/moy

    return error

def fast_neigh(data,tup,d_max=None,b=1,rm_center=False):
    """
    Return the neighbourhood of a given data point for an arbitrary size numpy array
    For now, accept only scalar data.
    The neighbourhood is chosen along the first d_max dimensions. If d_max is not specified, it returns the neigh along
    all the dimensions of the array
    INPUT
    -----

    OUTPUT
    -----
    data_neigh : 1d list
        List of values of the neighboors.
        The order is then lost in this operation.
    """
    dim = data.shape
    d_neigh = len(tup)

    if d_max is None:
        d_avg = d_neigh
    else:
        d_avg = d_max

#    if len(dim)==d_neigh:
#        print('scalar quantity')
#    else:
#        print('vector quantity')
    neigh = tuple()
    for k in range(d_avg):
        index_min = max(tup[k]-b,0)
        index_max = min(tup[k]+b+1,dim[k])
        neigh += (slice(index_min,index_max),)

    #if d_neigh>d_avg:
    for k in range(d_avg,d_neigh):
        neigh += (tup[k],)

    data_neigh = data[neigh]
    dim_neigh = data_neigh.shape

   # print(dim_neigh)

    tup_shape =(np.prod(dim_neigh[0:d_neigh]),)
    tup_shape += tuple(dim_neigh[d_neigh:]) # add other dimensions
    data_neigh = np.reshape(data_neigh,tup_shape)

    data_neigh = data_neigh[~np.isnan(data_neigh)]
    data_neigh = np.ndarray.tolist(data_neigh)

    if rm_center:
        if not np.isnan(data[tup]): #otherwise, is has already be removed !!
            data_neigh.remove(data[tup])
#        except:
#            print(data[tup])
    return data_neigh


def smoothn(dlist,s):
    #dlist is a list of numpy array to process
    #compute an interpolate version of d. Nan values are removed.
    #assume that z values are regularly spaced.
    # z is an n-dimension array ?
    #larger value of s means stronger smoothing
    di=[]

    for d in dlist:
        #Method 1 : use interpolate.Rbf
        dimensions=d.shape
        coord=np.where(~np.isnan(d))
        coord_i=np.where(d)

        #select automatically the dimension ??
        if len(dimensions)==1:
            rbfi = scipy.interpolate.Rbf(coord[0],d,function='multiquadric',smooth=s)
            di.append(rbfi(coord_i[0]))

        if len(dimensions)==2:
            rbfi = scipy.interpolate.Rbf(coord[0],coord[1],d,function='multiquadric',smooth=s)
            di.append(rbfi(coord_i[0],coord_i[1]))

    return di


def smooth(V,Dt):
    """
    average the data over a window Dt along its last axis
        The average is performed along the last axis of an N dimensionnal array (hopefully time axis)
    INPUT
    -----
    V : N dimensionnal numpy array
        The smoothing will be performed by default along the last axis
    Dt : int
        window width for the averaging

    OUTPUT
    -----
    Vs : numnpy array
        last axis has a shortened dimension by Dt
    """

    if Dt>1:
        dim=np.shape(V)
        nt=dim[-1]
        N=len(dim)

        Vs=np.array([np.nanmean(V[...,k-Dt:k+Dt],N-1) for k in range(Dt,nt-Dt)])
        Vs = np.transpose(Vs,tuple(range(1,N))+(0,))
    #some elements are lost during the smoothing operation : we should add the initial and enn elements
        Vs = np.concatenate((V[...,:Dt],Vs,V[...,-Dt:]),axis=N-1)
 #       Vs=np.transpose(np.array([V[...,:Dt]+[np.nanmean(V[...,k-Dt:k+Dt],N-1) for k in range(Dt,nt-Dt)]),axis=tuple(range(1,N))+(0,))
#    tup=
        print(Vs.shape)
        return Vs
    else:
        return V

### gaussian filter
def gaussian_smooth(M,field,sigma=0.5):
    data = getattr(M,field)
    data = filters.gaussian_filter(data, sigma=sigma, order=0, output=None)
    setattr(M,field,data)
    return M


def example_1():
    a=np.random.rand(10,10)
    t=(0,3)
    neighboors(a,t)

def example_2():
    n=100
    N=50
    sigma=N/5

    i=round(N/2)
    j=round(N/2)
    X=np.asarray([math.exp(-(k-i)**2/(2*sigma**2)) for k in range(N)])
    Y=np.asarray([math.exp(-(k-j)**2/(2*sigma**2)) for k in range(N)])

    Y=np.reshape(X,(N,1))
    X=np.reshape(X,(1,N))
    x=np.arange(N)
    y=np.arange(N)
    #generate
    a_init=np.dot(Y,X)

    i=np.floor(np.random.rand(n)*N)
    j=np.floor(np.random.rand(n)*N)

    a_dump=np.dot(Y,X)
    for t in zip(i,j):
        a_dump[t]=np.nan
    graphes.color_plot(x,y,a_dump,fignum=1)

    #remove NaN values
    rm_nans([a_dump])

    graphes.color_plot(x,y,a_dump,fignum=2)
    graphes.color_plot(x,y,a_dump-a_init,fignum=3)
    print(np.max(a_dump-a_init))

def examples():

    example_1()

    example_2()
    #second example : centered gaussian function
def trash():
    """
    delete some aberrant measurements before performing the iterations
    INPUT
    -----
    U : np array
    epsilon : float
        ??
    thresh : float
        ??
    OUTPUT
    -----
    U : np array
    """
    #parameters used for the filtering
    epsilon=0.02; #smaller than the noise value
    ny,nx=U.shape

    print(U.shape)
    nt=1

    #medianres=zeros(J,I);
    normfluct=np.zeros((ny,nx,2))
    med = np.zeros((ny,nx,nt))
    b=1;
    #eps=0.1;
    for t in range(nt):
        for i in range(nx):
            for j in range(ny):
                tup = (j,i,t)

                imin = max(i-b,0)
                imax = min(i+b+1,nx)
                jmin = max(j-b,0)
                jmax = min(j+b+1,ny)

                neigh = [(k,l) for l in range(imin,imax) for k in range(jmin,jmax)]
                neigh.remove()    #remove the center point
                U_neigh = [U[k] for k in neigh]
                print(neigh)
                print(U_neigh)
                #???    CHANGE THE FOLLOWING 5 Lines
#                neighcol=np.reshape(neigh,(n**2,)) # 8 neighbours + the central point
#                neighcol2=np.concatenate(neighcol[:n*b+b],neighcol[n*b+b+2:]) # list of the 8 neighbouhrs elements

                #compute the median of neighbouring elements
                med[j][i][t]=np.median(U_neigh)
                fluct=V[j][i]-med[j][i][t]
                res=neigh-med[j][i][t]
                medianres=np.median(abs(res));
                normfluct[j][i][t]=abs(fluct/(medianres+epsilon))

    #use logical index to set at NaN values the measurements exceeding a certain threshold
    #info1=(np.sqrt(normfluct[:][:][0]**2+normfluct[:][:][1]**2)>thresh)
#    utable[info1==1]=NaN;
 #   vtable[info1==1]=NaN;
    condition = (np.sqrt(np.sum(normfluct**2,2))>thresh)
    U = np.where(condition,med,U)

    return U
    #find typevector...
   # maskedpoints=numel(find((typevector)==0))
   # amountnans=numel(find(isnan(utable)==1))-maskedpoints;
    #discarded=amountnans/(size(utable,1)*size(utable,2))*100;
    #disp(['Discarded: ' num2str(amountnans) ' vectors = ' num2str(discarded) ' %'])


#main()
