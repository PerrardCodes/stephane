


"""
Store a dictionary into a hdf5 file
Look recursively for any dictionary contained
"""

import inspect
import numpy as np
import h5py
import os.path

import copy

#import stephane.manager.file_architecture as file_architecture

def convert(obj,L=[],t=0):
    """
    Convert recursively an homemade class to a dictionnary. 
    If a loop between homemade class instance is detected, cut the loop.
    INPUT
    -----
    
    OUTPUT
    -----
    """
    types = ['stephane','instance']
    
   # print(type(obj))
    if type(obj)==dict:
        return obj
    elif True in [a in str(type(obj)) for a in types]:
        name = copy.deepcopy(getattr(obj,'__name__')())
        if not name in L:
            L.append(name)            
            D = copy.deepcopy(getattr(obj,'__dict__'))
    
            for key in D.keys():
             #   print(key)
                if t<10:
                    attr = copy.deepcopy(getattr(obj,key))
                    D[key]=convert(attr,L=L,t=t+1)
                
            return D
        else:
            return 'pointer'
    else:
        return obj    
        
        

def load(obj,data):
    """
    Load a h5py group into the attributes of obj. If the h5py group contains subgroup, recursively saved them
    """
    if isinstance(data,h5py._hl.group.Group):        
        for key,item in data.attrs.items():
         #   print(obj)
            if type(obj)==dict:
                obj[key]=item
            else:
                setattr(obj,key,item)
                
        for key,item in data.items():
            if type(item) is h5py._hl.dataset.Dataset:
              #  print(key)
             #   print(obj)
                if type(obj)==dict:
                    obj[key]=item.value
                else:
                    setattr(obj,key,item.value)
                    
            elif type(item) is h5py._hl.group.Group:
                if not hasattr(obj,key):
                    setattr(obj,key,{})
                load(getattr(obj,key),item)

    else:
        print("Data are not a h5py group")
    return obj

############## Generating and Writing ####################   

def create(filename):
    import stephane.manager.file_architecture as file_architecture
    filename = file_architecture.os_i(filename)    
    
    Dirname = os.path.dirname(filename)
    if not os.path.exists(Dirname):
        os.makedirs(Dirname)        
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
    dict_total = convert(obj,L=[])
#    dict_total = get_attr_rec({},obj,[])

    if filename is None:
        filename = os.path.dirname(obj.Sdata.fileCine) + '/hdf5/test'+'.hdf5'         
    
    f,do = create(filename)
    if do:
        write_rec_terminal(f,dict_total,L=[],key=key)
        f.close()
        

def write_rec(f,Dict,key='',group=None,t=0,tmax=50):
    """
    Write recursively a dictionnary into a h5py previously open file (f)
    """
    
    done = False
    if type(Dict)==dict:
        done=True
        if group is None:
            group = key
        else:
            group = group+'/'+key
        
        if group not in f:
            f.create_group(group)
        for key in Dict.keys():
            write_rec(f,Dict[key],key=key,group=group)

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
            
    if type(Dict) in [int,str,float,np.int64,np.float64,bool]:
        done=True
        f[group].attrs[key] = Dict
    if not done:
        print("Unrecognized : "+str(key) +' of type '+str(type(Dict))) 
        
def write_rec_terminal(f,Dict,L=[],key='',group=None,t=0):
    """
    Write recursively a dictionnary into a h5py previously open file (f)
    """
    
    done = False
    if type(Dict)==dict:
        done=True

        if Dict in L:
            print('Already written ! recursive loop detected')
        else:
            L.append(Dict)
            if group is None:
                group = key
            else:
                group = group+'/'+key
        
            if group not in f:
               # print(group)
                f.create_group(group)
            for key in Dict.keys():
                write_rec_terminal(f,Dict[key],L=L,key=key,group=group,t=t+1)

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
            
    if type(Dict) in [int,str,float,np.int64,np.float64,bool]:
        done=True
       # print(key)
        f[group].attrs[key] = Dict
    if not done:
        print("Unrecognized : "+str(key) +' of type '+str(type(Dict))) 
        
     
############### Open and load ###########


def open(filename,typ='r'):
    """
    Open a hdf5 file. 
    """
   # filename = file_architecture.os_i(filename)            
    if os.path.exists(filename):
        f = h5py.File(filename,typ)
    else:
        print("File "+filename+" does not exist")
        f = h5py.File(filename,'w')
    return f
   # print('done')

def load_dict(data):
    """
    Transform a h5py group into a dictionnary. Recursive function
    """
    return load_rec(data) 

def load_rec(data):
    ans = {}
    for key,item in data.attrs.items():
        ans[key] = item         
        
    for key,item in data.items():
        if type(item) is h5py._hl.dataset.Dataset:
            ans[key] = item.value
        elif type(item) is h5py._hl.group.Group:
            ans[key] = load_dict(item)
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

def example(filename):    
    f,b = create(filename)
    
    Dict = {'Data':np.zeros(10),'param':{'p':3,'toto':True}}
    if b:
        write_rec_terminal(f,Dict,key='U')
        
def example_2(filename,Dict,key='U'):    
#    Dict = {'Data':np.zeros(10),'param':{'p':3,'toto':True,'subparameter':{'a':0,'b':np.random.random(20)}}}
    f,b = create(filename)   #create the file. b=True if it worked
    
    if b:
        write_rec_terminal(f,Dict,L=[],key=key)  #write recursively Dict into f
        
    f = open(filename)   # open an hdf5 file
    
    data = load_dict(f)  # load recursively f into data.
    #data has the same shape as the initial Dict
    
    return data  