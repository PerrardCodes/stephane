


import os
import numpy as np
import glob
import os.path



num = 170
folder='/Volumes/Stephane_Data_1/Stephane_Data_1/Experiments/2016_11_30/PIV_data/'#PIVlab_ratio2_W32pix_Dt_1_PIV_sv_vp_X10mm_Makro100mm_tube_fps5000_A0mm_piston12mm_v250mm_f5Hz__Cam_12621_Cine5/'

folders = np.asarray(glob.glob(folder+'*'))


for folder in folders:
    l= glob.glob(folder+'/*1.txt')
#    print("")
#    print("")
#    print("")
    print(l)
#    print(folder)
    #remove_short_name(l)
    for name in l:
        info = os.stat(name)
        if info.st_size<100000 or info.st_size>200000:
            print(info.st_size)
        os.remove(name)

def remove_short_name(l):
    l = np.asarray(l)
    u = np.asarray([len(i) for i in l])

    print(min(u))
    print(max(u))

    num = min(u)
    rm_files = l[np.where(u==num)]

    if max(u)>=num+2:
        print(rm_files)
        for name in rm_files:
            os.remove(name)

    
