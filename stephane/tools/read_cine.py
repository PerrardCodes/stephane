# -*- coding: utf-8 -*-
"""
Created on Tue Mar 17 15:38:03 2015

@author: stephane
"""
import stephane.Image_processing.cine as cine
import numpy as np
import pylab as plt
import os
import os.path

def read_cine(fileCine,n,debut=0):
    #n :number of frame to be read
    print('cine open')
    c = cine.Cine(fileCine)

    frameList=[]
#    print('toto')
    for i in range(debut,n):
        frame=read_frame(c,i)        
        frameList.append(frame)
        
    c.close()
    print('cine closed')
    
    return frameList   
    
def read_frame(cine,i):
        # get maximum value according to bit depth
        bitmax = float(2**cine.real_bpp)
        # get frames from cine file
        a = cine[i].astype("f")
        a = a / bitmax * 255
        
        im_ref = np.float64(a)
        
        return im_ref
        
def save_fig(fignumber,filename,Dir='',fileFormat='pdf'):
    if not Dir=='':    
        if not os.path.isdir(Dir):
            os.makedirs(Dir)
    # define file name
#        fileName = Dir + "n" + str(num) + '.' + fileFormat
    filename=filename+'.'+fileFormat
#    print(filename)
    # make the right figure active
    plt.figure(fignumber)    
    # save the figure
    dpi=50
    plt.savefig(filename, dpi=dpi)