import matplotlib.pyplot as plt
import glob
import numpy as np
import fluids2d.piv as piv
import stephane.cine.cine as cine
import sys
from functools import partial
from multiprocessing import Process,Pool
import os

#!/usr/bin/python
# -*- coding: utf8 -*-

from Mesure import Mesure
import pre_traitement as pre

import fluids2d.piv as p
#from danjruth.piv import *
from functools import partial
from multiprocessing import Process, Pool
import os
import numpy as np

class Volume(Mesure):
    def __init__(self, data, m={}):
        Mesure.__init__(self, data)
        self.m=m

    def analysis(self, parent_folder, cine_name, adresse_s, npy=None, fx=1., dt_origin="", \
                frame_diff="", crop_lims=None, maskers=None,\
                window_size=32, overlap=16, search_area_size=32,\
                save=True, s2n_thresh=1.2, bg_n_frames=None, a_frames=""):

        #Crée un objet processing présent dans : danjruth.piv
        processing = p.PIVDataProcessing(parent_folder, cine_name, name_for_save=adresse_s, dx=fx, dt_orig=dt_origin,\
                                        frame_diff=frame_diff, crop_lims=crop_lims, maskers=maskers,\
                                        window_size=window_size, overlap=16, search_area_size=search_area_size)
        #ajoute à m ce qu'il a dans l'objet processing
        self.m.update(processing.__dict__)
        #Si npy n'est pas donné en paramètre il lance l'analyse sur l'objet
        if(npy==None):
            flowfield = processing.run_analysis(a_frames=a_frames, save=save, s2n_thresh=s2n_thresh, bg_n_frames=bg_n_frames)
        #Sinon il load juste le npy dans m
        else :
            self.m['U'] = np.load(npy)

        return flowfield


    def load(self, npy):
        self.m['U'] = np.load(npy)

    def analysis_multi_proc(self, parent_folder, cine_name, adresse_s, npy=None, fx=1., dt_origin="", frame_diff="", crop_lims=None, maskers=None, window_size=32, overlap=16, search_area_size=32,save=True, s2n_thresh=1.2, bg_n_frames=None):
        Ncpu = os.cpu_count()
        with Pool(processes=Ncpu) as pool:
            ite = []
            for i in range(frame_diff):
                ite.append(np.arange(i, 127000-frame_diff, frame_diff))
	        #ite = [(0,25), (26, 50), (51, 75), (76, 100)]
            func = partial(self.analysis, parent_folder, cine_name, adresse_s, npy, fx, dt_origin,frame_diff, crop_lims, maskers, window_size, overlap, search_area_size, save, s2n_thresh, bg_n_frames)
            f = pool.map(func, ite)

            #get the image shape and the number of images processed by cpu from the output
            N=0
            for i in range(Ncpu):
                dim = f[i].shape
                N += dim[0]
            dimensions = (N,)+dim[1:] #assume all the images have the same shape
            flowfield = np.zeros(dimensions)

            for i in range(Ncpu):
                flowfield[i:N:Ncpu,...]=f[i]
            np.save(parent_folder+adresse_s+'_flowfield.npy',flowfield)
        self.m['U'] = flowfield
        return self


    def add_measurement(self, obj, name):
        setattr(self, obj, name)

    def get_name(self):
        return "PIV3D"



def find_timejumps(c,dtmin,dtmax):
    N = len(c)

    t =   np.asarray([c.get_time(i) for i in range(N)])
    jumps = np.logical_and(np.diff(t)>dtmin,np.diff(t)<dtmax)
    indices = np.where(jumps)[0]
    waittime = np.diff(t[indices])
    print('Maximum difference between waiting times : '+str(np.max(waittime)-np.min(waittime)))

    instant = []
    for i,ind in enumerate(indices[1:-1]):
        start = indices[i]+1
        end = indices[i+1]-1
        #plt.plot(np.diff(timages[start:end]))
        instant.append((start,end))
    return instant


def analysis_multi_proc(c,instant):
    #cut_function
    Ncpu = os.cpu_count()
    ite = []#np.zeros(Ncpu)
    instant = np.asarray(instant)
    for i in range(Ncpu):
        ite.append(np.arange(i,N,Ncpu))
    print(len(ite))
    print(ite[0])

    with Pool(processes=Ncpu) as pool:
	        #ite = [(0,25), (26, 50), (51, 75), (76, 100)]
            func = partial(find_positionlaser,c)
            f = pool.map(func, ite)

            instantV = []
            tV = []
            for i in range(Ncpu):
                instantV = instantV + f[i][0]
                tV = tV + f[i][1]

    return (instantV,tV)

def find_positionlaser(c,instant,a=5,Nz=None):
    instantV = []
    tV = []

    for (start,end) in instant:
        C = []
        print(start)
        for i in range(start,end):
            mean1 = np.mean(c.get_frame(i),axis=(0,1))
            mean2 = np.mean(c.get_frame(i+1),axis=(0,1))
            std1 = np.std(c.get_frame(i),axis=(0,1))
            std2 = np.std(c.get_frame(i+1),axis=(0,1))

            C.append(np.mean((c.get_frame(i)-mean1)*(c.get_frame(i+1)-mean2),axis=(0,1))/(std1*std2))

        maximum=[]
        minimum=[]
        for i in range(a,len(C)-a):
            window = slice(i-a,i+a+1)
            if np.argmax(C[window])==a+1:
                maximum.append(i+1)
                #plt.plot(i+1-maximum[0],C[i+1],'rx')
                if len(maximum)>1 and len(minimum)>0:
                    # get which way we are scanning
                    if (minimum[-1]-maximum[-2])<=Nz/2:
                        startV = maximum[-2]+start
                        endV = maximum[-1]+start
                    else:
                        startV = maximum[-1]+start
                        endV = maximum[-2]+start
                    instantV.append((startV,endV))
                    tV.append(c.get_time((startV+endV)//2))

            if np.argmin(C[window])==a+1:
                if len(maximum)>0:
                    minimum.append(i+1)
                    #if (minimum[-1]-maximum[-1])<=Nz/2:
                    #    plt.plot(i+1-maximum[0],C[i+1],'bo')
                    #else:
                    #    plt.plot(i+1-maximum[0],C[i+1],'k*')
                        #do a mirror symetry

        #plt.plot(C[maximum[0]:maximum[-1]])
        #plt.axis([0,150,0.65,1])
    return (instantV,tV)


if sys.platform=='win32':
    base = 'F:'
if sys.platform=='linux':
    base = '/media/stephane/DATA'
if sys.platform=='darwin':
    base = '/Volumes'

date = '20181106'
folder = base+'/Experimental_data/Turbulence3d/'+date+'/'
l =glob.glob(folder+'*.cine')
for i,name in enumerate(l):
    print(str(i)+' : '+os.path.basename(name))

#s = input()
s=3
try:
    i = int(s)
except:
    print("cannot be converted to an integer")

cinefile = l[i]
c = cine.Cine(cinefile)

#detecte les débuts et fin de Volumes
ft = 1./40000 #should be given directly by Data.param.ft
dtmin = 10*ft #we look for jumps at least 10 times Dt
dtmax = 10 #value in second
instant = find_timejumps(c,dtmin,dtmax)
N = len(c)

Nz = 25  #we should find the number of images using framerate/f, see data.param, need lea's package
#(instantV,tV) = find_positionlaser(c,instant,Nz=Nz)
(instantV,tV) = analysis_multi_proc(c,instant)
