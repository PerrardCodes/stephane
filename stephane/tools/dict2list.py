
import numpy as np

def add_list_to_dict(a_t,a):
    #concatene a primary dictionnaray (a_t) with a new one (a) in which each field of a_t (and a) are lists
    for key in a.keys():
        if key in a_t.keys():
            a_t[key]+=a[key]
        else:
          #  print("No "+key+" key, create a new dict key")
            a_t[key]=a[key]
    return a_t
    
def dict_to_dictlist(d):
    for key in d.keys():
        d[key] = to_1d_list(d[key])
    return d

def to_1d_list(a,d=3,t=False):
    dimensions = a.shape    
    if len(dimensions)==d:  
        return np.ndarray.tolist(np.reshape(a,(np.prod(a.shape),)))
    else:
        if len(dimensions)==1:
            return np.ndarray.tolist(a)
        else:
            return np.ndarray.tolist(np.reshape(a,(np.prod(a.shape[:d]),a.shape[-1])))           
         

def extract_dict(d,indices):
    d_extract = {}
    for key in d.keys():
        d_extract[key] = np.asarray(d[key])[indices]
    return d_extract