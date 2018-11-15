



import time
import copy



def make(M,field,step):
    """
    Use to corse grain the velocity fields
    """
M_corse = []
steps = range(1,50)
t1=time.time()
for step in steps:
    print(step)
    M = get_jhtd.JHTDdata(data,param,N=128)
    vgradient.compute(M,'omega',type=2,step=step)
    setattr(M,'step',step)
    M_corse.append(M)
t2=time.time()
print("Elapsed time : "+str(t2-t1))


def make_list(M,field,steps,**kwargs):
    key = field+'_grain'
    D = {}
    setattr(M,key,D)
    
    for step in steps:        
        data = measure(M,field,step=step,**kwargs)[0]
        D[data] = 