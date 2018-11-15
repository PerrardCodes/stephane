


"""
Store an instance of a homemade class to a json list of dictionnary. 
Look recursively for any attributes that are homemade classes.
Stop the loop when it goes back to a previously stored object (in case of reciprocal link between objects)
"""

import inspect
import json
import numpy as np

import os.path

import stephane.manager.file_architecture as file_architecture

def get_attributes(obj):
    dict_attr = {}
    recursive = []
    List = dir(obj)
    
    for arg in List:
        unfound = True
        if not arg[0:2]=='__':
            elem = getattr(obj,arg)
            if isinstance(elem,type(obj)):
                unfound = False
                recursive.append(elem)
             #   print('recursive, '+str(arg) + " : "+ str(type(elem)))
            if inspect.ismethod(elem):
                unfound = False
          #  print('method, '+str(arg) + " : "+ str(type(elem)))
            if unfound:
                if type(elem) in [dict,str,float,int,np.int64,np.float64]:
                #    print(type(elem))
         #   print('attribute, '+str(arg) + " : "+ str(type(elem)))
                    dict_attr.update({str(arg):elem})
        else:
            #not a generated attribute 
            pass
     
    return dict_attr,recursive
    
def get_attr_rec(dict_obj,obj,obj_list,convert=True):
    if not obj in obj_list:
        obj_list.append(obj)
        dict_present,recursive = get_attributes(obj)
        dict_obj[str(obj)] = dict_present
        
        for rec_obj in recursive:
         #   print(rec_obj)
            dict_par = get_attr_rec(dict_obj,rec_obj,obj_list)
            if not dict_par is None:
                dict_obj.update(dict_par)
        
        return dict_obj
    else:
        return None
            
def write_attributes(obj):
    dict_attr = get_attributes(obj)[0]
#    print(dict_attr)
    dict_attr = convert_array(dict_attr)
    #print(dict_attr)
    
def convert_array(dict_attr,erase=False):
    for key in dict_attr.keys():
        #print(dict_attr[key])
        if type(dict_attr[key])==np.ndarray:
            if erase==False:
            #    print(key+" converted from np array to List")
                dict_attr[key] = np.ndarray.tolist(dict_attr[key])
            else: 
                dict_attr[key] = 'Data_removed'
        if type(dict_attr[key])==dict:
            dict_attr[key] = convert_array(dict_attr[key],erase=erase)
       # print(type(dict_attr[key]))        
    return dict_attr
    
    
def write_rec(obj,erase=False,filename=None):
    """
    Write into a json file all the parameters recursively
    json file contains a dictionnary for each Class instance (e.g. Sdata, param, id) the parameters
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
#    for key in dict_total.keys():
#        print(key,dict_total[key])       
    for key in dict_total.keys():
        dict_total[key] = convert_array(dict_total[key],erase=erase)
        
    json_format = json.dumps(dict_total)

        
    if filename is None:
        filename = obj.param.fileParam[:-4] + '.json'
    
    Dir = file_architecture.os_c(os.path.dirname(filename))
    if not os.path.isdir(Dir):
        os.makedirs(Dir)
    filename = Dir + '/' +os.path.basename(filename)
    
    print(filename)
    f = open(filename,'w')
    f.write(json_format)
    f.close()
    
    print('done')
    
    
    