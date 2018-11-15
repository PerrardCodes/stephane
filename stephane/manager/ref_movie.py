# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 15:27:39 2015

@author: stephane
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Mar 12 14:42:10 2015

@author: stephane
"""

import glob
import os.path

import stephane.manager.cine2pic as cine2pic
import stephane.tools.browse as browse
import stephane.tools.rw_data as rw_data
#generate an small avi movie for each of the cine file, only the first 200 images to be used as spatial reference


def gen(directory,nimage=200):
    cineList = glob.glob(directory+'*.cine')
    #make a dir for references movies
    make(cineList,nimage=nimage)
        
def make(fileList,start=0,stop=200,nimage=20):       
    for file in fileList:
        folder = os.path.dirname(file) + '/References'
        make_movie(file,nimage=nimage,start=start,stop=stop,folder=folder)
    
def make_ref(fileList):
    for file in fileList:
        folder = os.path.dirname(file) + '/References'
        make_movie(file,nimage=20,start=0,stop=200,folder=folder)
        
def make_movie(cinefile,nimage,start=0,stop=1,step=1,folder='',framerate=30,quality=50):
    #generate a folder with tiff images : it would be easier if a avi file is generated directly
    root = os.path.basename(cinefile)
    base = browse.get_string(root,'',end='.cine')
    
    if folder=='':
        folder = os.path.dirname(cinefile)

    folder = '/References'
    print(folder+base)

    cine2pic.cine2tiff(cinefile,mode='Sample',step=step,ctime=2,start=start,stop=nimage*step+start,folder=folder,post=base)
    
#    indexList,ndigit = cine2pic.manual_sampling([0],[nimage])    
    
def make_ref_file(cinefile,folder=''):
    keys = ['fx','im0','x0','y0','angle']
    List_info = [[0] for i in range(len(keys))]
    
    if folder=='':
        folder = os.path.dirname(cinefile)

    #generate an empty ref file associated to each cine file 
    #Should be manually filled out from images measurements
    
    cinebase = browse.get_string(os.path.basename(cinefile),'','.cine')
    name = folder +'Ref_'+cinebase+'.txt'
    
    if not os.path.isfile(name):    
        print(name)
    #check first if a Ref file exists already    
#    Ref_PIV_sv_vp_zoom_Polymer_200ppm_X25mm_fps5000_n18000_beta500mu_H1180mm_S300mm
        rw_data.write_dictionnary(name,keys,List_info)
    else:
        print('Reference file exists already')
        pass
    
#directory = '/Volumes/labshared3/Stephane/Experiments/Accelerated_grid/2015_03_24/'
#gen(directory)