

import copy

def convert(obj,L=[],t=0):
    types = ['stephane','instance']
    
    print(type(obj))
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
        
        
def convert_manual(obj):
    #Sdata :
    D = copy.deepcopy(obj.__dict__)
    D['Id']=copy.deepcopy(obj.Id.__dict__)
    D['Id']['Sdata']=D

    return D