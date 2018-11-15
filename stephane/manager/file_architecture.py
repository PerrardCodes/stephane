# -*- coding: utf-8 -*-
"""
Created on Thu Oct  1 11:54:32 2015
@author: stephane
This file list all the possible locations of my data (Stephane, Gautam, Some of Joseph, examples from Martin) 
If a new directory tree is created for storing data, just add its root on the list below (Dir_data variable)
If you get trouble to access the data from a distant computer, 
add a new extension "distant_root" to the list below and an additionnal case in the method os_c
"""

#return the folders associated to each date, and other manipulations of files from different location
import os.path
import sys
import numpy as np

Dir_data=[]

distant_root = '/Volumes/'

Dir_data.append(distant_root+'labshared3/Stephane/Experiments/Accelerated_grid/')
Dir_data.append(distant_root+'labshared/Stephane_lab1/')
Dir_data.append(distant_root+'labshared/Stephane_lab1/Vortices/')
Dir_data.append(distant_root+'labshared/Stephane_lab1/Vortex_vs_turbulence/')


Dir_data.append(distant_root+'labshared2/Bipul/')
#Dir_data.append(distant_root+'labshared2/Weerapat/Four_vortex_collider/')
			
#Dir_data.append(distant_root+'Stephane/Vortices/Stephane/')

#Dir_data.append(distant_root+'labshared/Stephane_lab1/Vortices/') -> to be analyzed later

Dir_data.append(distant_root+'labshared/Gautam/Experiments/')

#Dir_data.append(distant_root+'labshared/Stephane_lab1/Vortex_vs_turbulence/') -> entirely moved to Irvinator-2

Dir_data.append(distant_root+'Stephane/Vortex_Turbulence/')
Dir_data.append(distant_root+'Stephane/Accelerated_grid/')
Dir_data.append(distant_root+'Stephane/Vortices/Gautam/')
Dir_data.append(distant_root+'Stephane/Vortices/Martin/')
Dir_data.append(distant_root+'Stephane/Vortices/Stephane/')

Dir_data.append(distant_root+'Stephane_Data_1/Stephane_Data_1/Experiments/')
#Dir_data.append(distant_root+'Stephane_Data_1-1/Experiments/')



local_root = '/Users/stephane/Documents/'

Dir_data.append(local_root+'Experiences_local/Accelerated_grid/')
Dir_data.append(local_root+'Experiences_local/Vortex_Turbulence/')
Dir_data.append(local_root+'Experiences_local/Others/')
Dir_data.append(local_root+'Experiences_local/Ballooon/')

Dir_data.append(local_root+'Experiences_local/Vortex_collider/')



Dir_data.append(local_root+'Gautam/')


linux_root = '/home/steph/Documents/'

Dir_data.append(linux_root+'Stephane/Vortex_Turbulence/')
Dir_data.append(linux_root+'Stephane/Accelerated_grid/')
Dir_data.append(linux_root+'Stephane/Vortices/Gautam/')
Dir_data.append(linux_root+'Stephane/Vortices/Martin/')

linux_root_2 = '/media/steph/'
Dir_data.append(linux_root_2+'Data_1/Stephane_Data_1/Experiments/')
#Dir_data.append(linux_root_2+'Stephane_Data_1-1/Experiments/')

#Dir_data.append(linux_root+'Stephane/Vortices/Stephane/')



def get_dir(date):
    """
    Return the directory in which the data corresponding to the date given are stored.
    Look in any of the specified rootfolder
    INPUT 
    -----
    date : str
        date to be processed
    OUTPUT
    path : str
        path of the rootdir
    -----
    """
    Dirdata=[]
    
    for Dir in Dir_data:
        if os.path.isdir(Dir+date):
            Dirdata.append(Dir+date)
    
    if Dirdata==[]:        
        print('No folder found for '+date)
        return ''
    if len(Dirdata)==1:
        return Dirdata[0]
    if len(Dirdata)>1:
#        print('Multiple choice of directory, make a choice among :')
       # for i,Dir in enumerate(Dirdata):
    #        print(str(i)+ ' : '+Dir)
#        j=0
        search = [distant_root+'Stephane' in elem for elem in Dirdata]
    #    print(search)
        if any(search):
            j = np.where(search)[0][0]
      #      print(j)
       #     input()
        else:
            j=0  #do not ask for input anymore. preset by the user, or take the first one of the list as a default one
#            j=int(input())
        return Dirdata[j]
        
        
def os_c(Dir):
   # print("Directory :" +Dir)
    
    if sys.platform=='darwin':
        if linux_root in Dir:
            n = len(linux_root)
            return distant_root+Dir[n:]
        if linux_root_2 in Dir:
            n = len(linux_root_2)
            print(distant_root+Dir[n:])
            return distant_root+Dir[n:]
            
    if (distant_root in Dir) and 'linux' in sys.platform:
        n = len(distant_root)
        return linux_root+Dir[n:]
    return Dir    
    
def os_i(filename):
    Dir = os_c(os.path.dirname(filename))
    if not os.path.isdir(Dir):
        os.makedirs(Dir)
    filename = Dir + '/' +os.path.basename(filename)
    return filename
        
        
        
        
        
        
        
        