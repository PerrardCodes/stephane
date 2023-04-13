########################################################################
#
#  Copyright 2014 Johns Hopkins University
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# Contact: turbulence@pha.jhu.edu
# Website: http://turbulence.pha.jhu.edu/
#
########################################################################

import os
import sys
import datetime
if sys.version_info[0] == 2:
    import urllib
elif sys.version_info[0] == 3:
    import urllib.request, urllib.parse, urllib.error

import time
import numpy as np
import h5py
import stephane.tools.rw_data as rw_data

def get_cutout(
        filename = 'test',
        t0 = 0, tl = 1,
        x0 = 0, xl = 16,
        y0 = 0, yl = 16,
        z0 = 0, zl = 16,
        data_set = 'isotropic1024coarse',
        data_type = 'u',
#        auth_token = 'edu.jhu.pha.turbulence.testing-201311',
        auth_token =  'edu.uchicago.sperrard-9cf78a64',
        base_website = 'dsp033.pha.jhu.edu'):
    url = ('http://{0}/jhtdb/getcutout/{1}/{2}/{3}/'.format(
                base_website, auth_token, data_set, data_type)
         + '{0},{1}/'.format(t0, tl)
         + '{0},{1}/'.format(x0, xl)
         + '{0},{1}/'.format(y0, yl)
         + '{0},{1}/'.format(z0, zl))
    """
# auth_token =  'edu.uchicago.sperrard-9cf78a64',
        base_website = 'turbulence.pha.jhu.edu'):
    url = ('http://{0}/cutout/download.aspx/{1}/{2}/{3}/'.format(
                base_website, auth_token, data_set, data_type)
"""
#http://dsp033.pha.jhu.edu/jhtdb/getcutout/
#edu.uchicago.sperrard-9cf78a64/isotropic1024coarse/u/0,16/0,16/0,16/16,1/hdf5/
         
    print(url)
    if data_type in ['u', 'b', 'a']:
        ncomponents = 3
    elif data_type in ['p']:
        ncomponents = 1
    elif data_type in ['ub']:
        ncomponents = 6
#    print('Retrieving h5 file, size {0} MB = {1} MiB.'.format(
#            xl*yl*zl*ncomponents * 4. / 10**6,
#            xl*yl*zl*ncomponents * 4. / 2**20))
    if os.path.isfile(filename+'.h5'):
        os.remove(filename+'.h5')
    
    if sys.version_info[0] == 2:
        urllib.urlretrieve(url, filename + '.hdf5')
    elif sys.version_info[0] == 3:
        urllib.request.urlretrieve(url, filename + '.hdf5')
    # check if file downloaded ok
    print(filename)
    data = h5py.File(filename + '.hdf5', mode = 'r')
    data.close()
    print('Data downloaded and ' + filename + '.h5 written successfuly.')
    return None

def get_big_cutout(
        filename = 'tst',
        t0 = 0, tl = 1,
        x0 = 0, xl = 32,
        y0 = 0, yl = 32,
        z0 = 0, zl = 32,
        chunk_xdim = 16,
        chunk_ydim = 16,
        chunk_zdim = 16,
        data_set = 'isotropic1024coarse',
        data_type = 'u',
        auth_token = 'edu.uchicago.sperrard-9cf78a64',
        base_website = 'turbulence.pha.jhu.edu'):
    big_data_file = h5py.File(filename + '.h5', mode='w')
    xchunk_list = [chunk_xdim for n in range(int(xl / chunk_xdim))]
    if not (xl % chunk_xdim == 0):
        xchunk_list.append(xl % chunk_xdim)
    ychunk_list = [chunk_ydim for n in range(int(yl / chunk_ydim))]
    if not (yl % chunk_ydim == 0):
        ychunk_list.append(yl % chunk_ydim)
    zchunk_list = [chunk_zdim for n in range(int(zl / chunk_zdim))]
    if not (zl % chunk_zdim == 0):
        zchunk_list.append(zl % chunk_zdim)
    for current_data_type in data_type:
        if current_data_type in ['u', 'b', 'a']:
            ncomponents = 3
        elif current_data_type in ['p']:
            ncomponents = 1
        big_data = []
        for time in range(t0, t0 + tl):
            big_data.append(big_data_file.create_dataset(
                    current_data_type + '{0:0>5}'.format(time*10),
                    (zl, yl, xl, ncomponents),
                    np.float32,
                    compression = 'lzf')) ### is compression a good idea?
        N = len(zchunk_list)*len(ychunk_list)*len(xchunk_list)*tl
        count=0
        for cz in range(len(zchunk_list)):
            for cy in range(len(ychunk_list)):
                for cx in range(len(xchunk_list)):
                    for time in range(t0, t0 + tl):
                        tmp_filename = (filename
                                      + '_{0:0>2x}{1:0>2x}{2:0>2x}_{3}'.format(cz, cy, cx, current_data_type))
                        
                        try:
                            count+=1
                            print(str(count*100/N)+' %')
                            if not os.path.exists(tmp_filename + '.h5'):
                                get_cutout(
                                        tmp_filename,
                                        t0 = time, tl = tl,
                                        x0 = x0+cx*chunk_xdim, y0 = y0+cy*chunk_ydim, z0 = z0+cz*chunk_zdim,
                                        xl = xchunk_list[cx], yl = ychunk_list[cy], zl = zchunk_list[cz],
                                        data_set = data_set,
                                        data_type = current_data_type,
                                        auth_token = auth_token,
                                        base_website = base_website)
                                    
                            new_file = h5py.File(tmp_filename + '.h5', mode='r')
                            new_data = new_file[current_data_type + '{0:0>5}'.format(time*10)]
                            big_data[time - t0][cz*chunk_zdim:cz*chunk_zdim+zchunk_list[cz],
                                            cy*chunk_ydim:cy*chunk_ydim+ychunk_list[cy],
                                            cx*chunk_xdim:cx*chunk_xdim+xchunk_list[cx], :] = new_data
                        except:
                            print('Data not loaded')
                            fail = np.empty((chunk_zdim,chunk_ydim,chunk_xdim,3))
                            fail.fill(np.nan)
                            big_data[time - t0][cz*chunk_zdim:cz*chunk_zdim+zchunk_list[cz],
                                            cy*chunk_ydim:cy*chunk_ydim+ychunk_list[cy],
                                            cx*chunk_xdim:cx*chunk_xdim+xchunk_list[cx], :] = fail
                        
    big_data_file.create_dataset(
            '_contents',
            new_file['_contents'].shape,
            new_file['_contents'].dtype)
    big_data_file['_contents'][:] = new_file['_contents'][:]
    if data_type == 'ub':
        big_data_file['_contents'][0] = 5
    big_data_file.create_dataset(
            '_dataset',
            new_file['_dataset'].shape,
            new_file['_dataset'].dtype)
    big_data_file['_dataset'][:] = new_file['_dataset'][:]
    big_data_file.create_dataset(
            '_size',
            new_file['_size'].shape,
            new_file['_size'].dtype)
    big_data_file['_size'][0] = new_file['_size'][0]
    big_data_file['_size'][1] = xl
    big_data_file['_size'][2] = yl
    big_data_file['_size'][3] = zl
    big_data_file.create_dataset(
            '_start',
            new_file['_start'].shape,
            new_file['_start'].dtype)
    big_data_file['_start'][0] = new_file['_start'][0]
    big_data_file['_start'][1] = x0
    big_data_file['_start'][2] = y0
    big_data_file['_start'][3] = z0
    new_file.close()
    big_data_file.close()
    return None
    
def convert_dict(dset,delimiter='_'):
    #convert a dictionnary to a string format, using 
    s=''
    for key in dset.keys():
        s+= key + delimiter + str(dset[key]) + delimiter

    return s[:-1]

def load_spatial_samples():
    # dumb test
    data_type = 'u'
    dset = dict(['t0','tl','x0','xl','y0','yl','z0','zl',])

    c0 = 50
    step = 100
    
    tstep = range(c0,1024,step)
    xstep = range(c0,1024,step)
    ystep = range(c0,1024,step)
    zstep = range(c0,1024,step)
    
    indices = [(t,i,j,k) for t in tstep for i in xstep for j in ystep for k in zstep]
        
    d=3
    N = len(tstep)**(d+1)
    
    nt = 1
    nx = 10
    ny = 10
    nz = 10
        
    print("Number of files to load : "+str(N))
    count = 0
    for i in range(N):
        print(i)
        keys = ['t0','x0','y0','z0','tl','xl','yl','zl']
        values = [indices[i][j] for j in range(d+1)]+[nt,nx,ny,nz]
        dset = {key : values[j] for j,key in enumerate(keys)}

        name = 'jhtdb_data_' + convert_dict(dset)# + '.h5'
        
        try:
            get_cutout(
                t0 = dset['t0'], tl = dset['tl'],
                x0 = dset['x0'], xl = dset['xl'],
                y0 = dset['y0'], yl = dset['yl'],
                z0 = dset['z0'], zl = dset['zl'],
                data_set = 'isotropic1024coarse',
                data_type = data_type,filename=name)
        except:
            count+=1
            pass
            
    print("Percentage of data lost : "+str(count*100//N)+ ' %')

      #  f0 = h5py.File(name+'.h5', mode='r')
      #  print((data_type, f0['_contents'][:]))
      #  f0.close()
    return None
    
def date():
    a = time.strptime(time.ctime())
    year = str(a[0])
    n1 = len(str(a[1]))
    n2 = len(str(a[2]))
    
    zeros = ''
    for i in range(2-n1):
        zeros+='0'
    month = zeros+str(a[1])
    zeros = ''
    for i in range(2-n2):
        zeros+='0'
    day = zeros+str(a[2])
    
    date = year+'_'+month+'_'+day
    return date
    
def main():
    dirbase = '/Users/stephane/Documents/JHT_Database/Data/Spatial_measurement_2d/'+date()+'/'#'#_2016_04_11/'   
    print(date())
#    input()
    
    directory = dirbase + 'Data/'
    keys = ['t0','x0','y0','z0','tl','xl','yl','zl']

    if not os.path.isdir(directory):
        print("Not a directory")
        os.makedirs(directory)
        
    N = 256    
    Nt = 1024
    N0 = 0
    t = range(Nt)#,2**4)#[256]#range(N/2**6)
    
    log_file = dirbase + 'log.txt';
    f_log = open(log_file,'w')
    begin = True
     
    for t0 in t:
        print(t0)
        values = [t0,512,0,0,1,1,N,N]#[indices[i][j] for j in range(d+1)]+[nt,nx,ny,nz]
        dset = {key : values[j] for j,key in enumerate(keys)}
        
        dir_current = directory + convert_dict(dset)+'/'
        if not os.path.isdir(dir_current):
            os.makedirs(dir_current)
        
        name = dir_current + 'jhtdb_data_' + convert_dict(dset)# + '.h5'
        data_type = 'u'
        print(name)
#    load_spatial_samples()
      #  try:
        get_cutout(
        t0 = dset['t0'], tl = dset['tl'],
        x0 = dset['x0'], xl = dset['xl'],
        y0 = dset['y0'], yl = dset['yl'],
        z0 = dset['z0'], zl = dset['zl'],
        data_set = 'isotropic1024coarse',
        data_type = data_type,filename=name)#,chunk_xdim = 1,chunk_ydim = 1,chunk_zdim = N/2)
        
        
        """
        except:
            if begin:
                begin = False
                List_key = dset.keys()
                rw_data.write_header(f_log,List_key)
            List = [dset[key] for key in List_key]
            rw_data.write_line(f_log,List)
            
            print('Failed')
            print(List)   
        """
    f_log.close()        
    
def recover(log_file):
    """
    From a log file containing the dset information to download, try again to download data
    """
    data_type = 'u'
    directory = os.path.dirname(log_file)+ 'Data/'
    print(directory)
    
    Header,dset_list = rw_data.read_dataFile(log_file,Hdelimiter='\t',Ddelimiter='\t',Oneline=False)
    
    n = len(dset_list[dset_list.keys()[0]])
    
    for i in range(n):
        dset = {key:int(dset_list[key][i]) for key in dset_list.keys()}
        print(dset)
        
        dir_current = directory + convert_dict(dset)+'/'
        if not os.path.isdir(dir_current):
            os.makedirs(dir_current)
        name = dir_current + 'jhtdb_data_' + convert_dict(dset)# + '.h5'
        
        get_cutout(
        t0 = dset['t0'], tl = dset['tl'],
        x0 = dset['x0'], xl = dset['xl'],
        y0 = dset['y0'], yl = dset['yl'],
        z0 = dset['z0'], zl = dset['zl'],
        data_set = 'isotropic1024coarse',
        data_type = data_type,filename=name)

    
if __name__ == '__main__':
#    log_file = '/Users/stephane/Documents/JHT_Database/Data/Spatial_measurement_1d_2016_03_30/log.txt'
    main()
#    recover(log_file)

