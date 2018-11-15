


"""
Store an instance of a homemade class to a json list of dictionnary. 
Look recursively for any attributes that are homemade classes.
Stop the loop when it goes back to a previously stored object (in case of reciprocal link between objects)
"""

import inspect
import numpy as np
import h5py
import os.path

import stephane.manager.file_architecture as file_architecture

def get_attributes(obj):
    f = Fake()
    
    dict_attr = {}
    recursive = []
    List = dir(obj)
    
    for arg in List:
     #   print(arg)
        unfound = True
      #  print(arg)
        if not arg[0:2]=='__':
            elem = getattr(obj,arg)
           # print(elem)
            
            if isinstance(elem,type(f)) or isinstance(elem,type(obj)) or 'Sdata' in str(type(elem)): 
                #Last boolean condition is bad, but how to do otherwise ???
                print(str(type(elem)))
                unfound = False
                recursive.append(elem)

            if inspect.ismethod(elem):
                unfound = False

            if unfound:
                dict_attr.update({str(arg):elem})
        else:
            pass
    recursive = []
    return dict_attr,recursive
    
def convert_dict():
    pass

#def get_attr(obj,name =['stephane','instance']):
    
#    if True in [a in str(type(obj)) for a in ['stephane','instance']]:
#        return 
#    else:
        

def get_attributes(obj):
    f = Fake()
    
    dict_attr = {}
    recursive = []
    List = dir(obj)
    
    for arg in List:
     #   print(arg)
        unfound = True
      #  print(arg)
        if not arg[0:2]=='__':
            elem = getattr(obj,arg)
           # print(elem)
            
            if isinstance(elem,type(f)) or isinstance(elem,type(obj)) or 'Sdata' in str(type(elem)): 
                #Last boolean condition is bad, but how to do otherwise ???
                print(str(type(elem)))
                unfound = False
                recursive.append(elem)

            if inspect.ismethod(elem):
                unfound = False

            if unfound:
                print(str(arg))
                dict_attr.update({str(arg):elem})
        else:
            pass
    recursive = []
    return dict_attr,recursive
    

    
def get_attr_rec(dict_obj,obj,obj_list,convert=True,t=0):
    
    if not obj in obj_list:
        obj_list.append(obj)
        dict_present,recursive = get_attributes(obj)
        
        dict_obj[str(obj)] = dict_present
        
        if t<1:
            for rec_obj in recursive:
         #   print(rec_obj)
                dict_par = get_attr_rec(dict_obj[str(obj)],rec_obj,obj_list,t=t+1)
                if not dict_par is None:
                    dict_obj[str(obj)].update(dict_par)
        
        return dict_obj
    else:
        return None
            
def write_attributes(obj):
    dict_attr = get_attributes(obj)[0]
   # print(dict_attr.__name__)
#    dict_attr = convert_array(dict_attr)
    #json.dumps(dict_attr)
    #print(dict_attr)
    
def convert_array(dict_attr,erase=False):
    for key in dict_attr.keys():
        if type(dict_attr[key])==type(np.zeros(0)):
            if erase==False:
            #    print(key+" converted from np array to List")
                dict_attr[key] = np.ndarray.tolist(dict_attr[key])
            else: 
                dict_attr[key] = 'Data_removed'
       # print(type(dict_attr[key]))        
    return dict_attr


def read(filename,class_obj):
    """
    Read a hdf5 and write it in the attributes given
    """
    f = h5py.File(filename,'r+')
    
    read_rec(f,class_obj)
    
def read_rec(f):
    
    if True:
        for key in f.keys():
            read_rec(f[key])
    else:
        pass
        
def load(obj,data):
    """
    Load a h5py group into the attributes of obj. If the h5py group contains subgroup, recursively saved them
    """
    if isinstance(data,h5py._hl.group.Group):
        for key,item in data.attrs.items():
            setattr(obj,key,item)
         #       print(key,getattr(self,key))
                
        for key,item in data.items():
            if type(item) is h5py._hl.dataset.Dataset:
                setattr(obj,key,item.value)
                #    print(key,getattr(self,key))
            elif type(item) is h5py._hl.group.Group:
                    #print(key,item)
               # print(key+" Contain subgroup ! iterate the generator")
                #    print(key)
                   # print(to_hdf5.load_dict(item))
                val = load_dict(item)
                #    print(val)
                setattr(obj,key,val)   # !!!!! For some reasons, the output is wrong. copy the value of the 
    else:
        print("Data are not a h5py group")
    return obj

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
    
def load_dict_from_hdf5(filename):
    """
    Load a dictionary whose contents are only strings, floats, ints,
    numpy arrays, and other dictionaries following this structure
    from an HDF5 file. These dictionaries can then be used to reconstruct
    ReportInterface subclass instances using the
    ReportInterface.__from_dict__() method.
    """
    with h5py.File(filename, 'r') as h5file:
        return recursively_load_dict_contents_from_group(h5file, '/')

def recursively_load_dict_contents_from_group(h5file, path):
    """
    Load contents of an HDF5 group. If further groups are encountered,
    treat them like dicts and continue to load them recursively.
    """
    ans = {}
    for key, item in h5file[path].items():
        if type(item) is h5py._hl.dataset.Dataset:
            ans[key] = item.value
        elif type(item) is h5py._hl.group.Group:
            ans[key] = recursively_load_dict_contents_from_group(h5file, path + key + '/')
        else:
            pass
           # print(item)
  #  for key, item in h5file[path].attrs.items():
       # print(key,item)
    return ans
    
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
        write_rec_terminal(f,dict_total,L=[],key=key)
        f.close()

def open(filename,typ='r'):
    """
    Open a hdf5 file. 
    Partially cross platform function : can switch the directory name between linux and mac syntax
    """
    filename = file_architecture.os_i(filename)            
    if os.path.exists(filename):
        f = h5py.File(filename,typ)
    else:
        print("File "+filename+" does not exist")
        f = h5py.File(filename,'w')
    return f
   # print('done')
   
def create(filename,overwrite=True):
    filename = file_architecture.os_i(filename)            
    if not os.path.exists(filename):
        f = h5py.File(filename,'w')
        return f,True
    elif overwrite:
        os.remove(filename)
        return create(filename)
    else:
        print("File "+filename+" already exists, skip ")
        return None,False
        
"""        
def write_rec(f,Dict,key='U'):
    
    if type(Dict)==dict:
        pass
    else:
        name = group+'/'+key
        write_data(f,Dict)
"""

def write_data(f,data,name):    
    if type(data) in [list]:
        Dict = np.asarray(data)
    
    if type(data) in [np.ndarray]:
         done=True
         if dataname not in f:
             dset = f.create_dataset(name, data = data, chunks = True)#Dict.shape, dtype=Dict.dtype)
         else:
             print('overwrite date')
             f[dataname][...] = data  #dimensions should already match !!!
        
    if type(data) in [int,str,float,np.int64,np.float64]:
        done=True
        # print(key)
        grp.attrs[key] = Dict
    if not done:
        print("Unrecognized : "+str(key) +' of type '+str(type(Dict)))
            
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
            print("f")
            print(key)
            print(f.keys())     
            print(group)
            s = group
            for key in ['instance','object']:
                if key in group:
                    s = group[:group.find(key)] #avoid duplication of instance elements stored at different memory locations.
            
      #      print(group,s)
            
            found = False
            for k in f:
                print("")
                print(k)
                print(s)
                if s in k:
                    found = True
                    grp = f[k]                
            if not found:
                grp = f.create_group(group)
        else:
            grp = f[group]
#        print(Dict.keys())            
        for key in Dict.keys():
            if t<tmax:#limit the number of recursion to 2 : protection against overflow
               # print(key)
                write_rec(f,Dict[key],key=key,group=group,grp=grp,t=t+1)

    if type(Dict) in [list]:
        Dict = np.asarray(Dict)
            
    if type(Dict) in [np.ndarray]:
        done=True
        dataname = group+'/'+key
        if dataname not in f:
            dset = f.create_dataset(dataname, data = Dict, chunks = True)#Dict.shape, dtype=Dict.dtype)
        else:
            f[dataname][...] = Dict  #dimensions should already match !!!
            
    if type(Dict) in [int,str,float,np.int64,np.float64]:
        done=True
       # print(key)
        grp.attrs[key] = Dict
    if not done:
        print("Unrecognized : "+str(key) +' of type '+str(type(Dict))) 
        
def display(Dict,key=None):
    #recursive display of a dictionnary
    if type(Dict)==dict:
        print("")
        print("iterate "+str(key))
        for key in Dict.keys():
            #print(key)
            display(Dict[key],key=key)
    else:
        print('     '+key+'  , ' + str(type(Dict)))

class Fake:
    def __init__(self):
        pass
        
def example():
    filename = './Test/test.hdf5'
    
    f,b = create(filename)
    
    Dict = {'Data':np.zeros(10),'param':{}}
    if b:
        write_rec(f,Dict)