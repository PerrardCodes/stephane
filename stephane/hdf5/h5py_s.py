


"""
Store a dictionary into a hdf5 file
Look recursively for any dictionary contained
"""

import inspect
import numpy as np
import h5py
import os.path

#import stephane.manager.file_architecture as file_architecture



############## Generating and Writing ####################   

def create(filename):
#    filename = file_architecture.os_i(filename)            
    if not os.path.exists(filename):
        print(filename)
        f = h5py.File(filename,'w')
        return f,True
    else:
        print("File "+filename+" already exists, skip ")
        return None,False
  
def write(obj,erase=False,filename=None,key=''):
    """
    Write into a hdf5 file all the parameters recursively
    hdf5 file contains a dictionnary for each Class instance (e.g. Sdata, param, id) the parameters
    each individual is a dictionnary containing the attributes of the class + '__module__' and '__doc__'
    INPUT
    -----
    obj : Class instance 
        to be writen in json file. Attribute can be any type, numpy array are replaced by List.
    erase : bool, default False
        erase numpy array data before writing the json file. Prevent large json files
    OUTPUT
    -----
    None
    
    """
    dict_total = get_attr_rec({},obj,[])

    if filename is None:
        filename = os.path.dirname(obj.Sdata.fileCine) + '/hdf5/test'+'.hdf5'         
    
    f,do = create(filename)
    if do: 
        write_rec(f,dict_total,key=key)
        f.close()
        

def write_rec(f,Dict,key='',grp=None,group=None,t=0,tmax=3):
    """
    Write recursively a dictionnary into a h5py previously open file (f)
    """
    done = False
    if type(Dict)==dict:
       # print("Write :" +str(Dict))
        done=True
        
        if group is None:
            group = key
        #    grp = f.create_group(group)
        else:
#            if 'object' in :
            group = group+'/'+key
        
        if group not in f:
            print(group)
            grp = f.create_group(group)
        else:
            grp = f[group]
#        print(Dict.keys())            
        for key in Dict.keys():
            if t<tmax:#limit the number of recursion to 2 : protection against overflow
                write_rec(f,Dict[key],key=key,group=group,grp=grp,t=t+1)

    if type(Dict) in [list]:
        done=True
        Dict = np.asarray(Dict)
            
    if type(Dict) in [np.ndarray]:
        done=True
        dataname = group+'/'+key
        if dataname not in f:
            dset = f.create_dataset(dataname, data = Dict, chunks = True)#Dict.shape, dtype=Dict.dtype)
        else:
            f[dataname][...] = Dict  #dimensions should already match !!!
            
    if type(Dict) in [bool,int,str,float,np.int64,np.float64]:
        done=True
       # print(key)
        grp.attrs[key] = Dict
    if not done:
        print("Unrecognized : "+str(key) +' of type '+str(type(Dict))) 
        
     
############### Open and load ###########



def open(filename,typ='r'):
    """
    Open a hdf5 file. 
    Partially cross platform function : can switch the directory name between linux and mac syntax
    """
   # filename = file_architecture.os_i(filename)            
    if os.path.exists(filename):
        f = h5py.File(filename,typ)
    else:
        print("File "+filename+" does not exist")
        f = h5py.File(filename,'w')
    return f
   # print('done')

def load_dict(data,ans={}):
    """
    Transform a h5py group into a dictionnary. Recursive function
    """
    for key,item in data.attrs.items():
        ans[key] = item         
    for key,item in data.items():
        if type(item) is h5py._hl.dataset.Dataset:
         #   print(item)
            ans[key] = item.value
                #    print(key,getattr(self,key))
        elif type(item) is h5py._hl.group.Group:
            ans[key] = load_dict(item,ans=ans)
#                    load(item,data,path='')
         #   print(key,item)
        #    print("Contain subgroup ! iterate the generator")
    return ans
    
      
        
def display(Dict,key=None):
    #recursive display of a dictionnary
    if type(Dict)==dict:
        for key in Dict.keys():
            display(Dict[key],key=key)
    else:
        print('     '+key+'  , ' + str(type(Dict)))

#################### How to use it ! #####################
        
def example(filename,Dict=None,key='U'):    
    f,b = create(filename)
    
    if Dict is None:
        Dict = {'Data':np.zeros(10),'param':0}
    if b:
        write_rec(f,Dict,key=key)
        
def example_2(filename,Dict,key='U'):
    f,b = create(filename)   #create the file. b=True if it worked
    
    if b:
        write_rec(f,Dict,key=key)  #write recursively Dict into f
        
    f = open(filename)   # open an hdf5 file
    
    data = load_dict(f)  # load recursively f into data.
    #data has the same shape as the initial Dict
    
    return data