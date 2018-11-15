# -*- coding: utf-8 -*-
"""
Created on Fri Sep  4 11:09:28 2015

@author: stephane
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Sep  4 10:39:19 2015

@author: stephane
"""

import numpy as np
import math

def average(X,axis=None):
    
    if axis == None:
        Xmoy=np.nanmedian(np.asarray(X))
        Xstd=np.nanstd(np.asarray(X))
    else:
        Xmoy=np.nanmedian(np.asarray(X),axis=axis)
        Xstd=np.nanstd(np.asarray(X),axis=axis)
        
    return Xmoy,Xstd
    
def box_average(tlist,N):
    # from a panel of curves, compute the average using an histogramm like technic (mean value of every data in a given range)
    #list of tuples, each of these containing (Xi,Yi) data (stored in list types ?) that will be used for the averaging
    # N : number of boxes used
    x=[]
    y=[]
    for tup in tlist:      
        x = x + np.ndarray.tolist(tup[0])
        y = y + np.ndarray.tolist(tup[1])

    X = np.asarray(x)
    Y = np.asarray(y)
    
    # sort the resulting arrays by x values (is it useful ?)
  #  indices = np.argsort(X)     
  #  X = X[indices]
  #  Y = Y[indices]     
   
    Dx = (max(x)-min(x))/N    
    xmin = np.arange(min(x),max(x)+Dx,Dx)[:-1]
    xmax = np.arange(min(x),max(x)+Dx,Dx)[1:]

    Xmoy = np.zeros(N)
    Xstd = np.zeros(N)
    Ymoy = np.zeros(N)
    Ystd = np.zeros(N)
    
    for i in range(N):
        #find all the values in a given range
        indices = np.logical_and(X>xmin[i],X<=xmax[i])    
  #      print("Number of elements :"+str(np.sum(indices)))
        
        Xmoy[i],Xstd[i]=average(X[indices])
        Ymoy[i],Ystd[i]=average(Y[indices])

    return Xmoy,Ymoy,Ystd        
        
        