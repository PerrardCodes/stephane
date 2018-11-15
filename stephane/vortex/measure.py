


"""
From a set of measurement (?), compute the following quantities :
- compute vorticity, strain and energy
- generate a movie for the first Mdata
- track the core. Show the position vs time along both axis
- Averaged distribution of vorticity over the entire movie. 
- Same for circulation

"""


import stephane.mdata.Sdata_manip as Sdata_manip 
import stephane.analysis.vgradient as vgradient

import stephane.display.graphes as graphes

date = '2016_11_12'
indices = range(60)
savedir = './Vortex_Turbulence/Vortex_propagation/'+date+'/'

def load(date,indices,mindices=0):
    Slist = Sdata_manip.load_serie(date,indices)
    Mlist = Sdata_manip.load_measures(Slist,indices=mindices)
    
    return Mlist

def compute(M):
    M.get('E')
    M.get('omega')
    M.get('strain')
    
def corse_grain(M,field='omega'):
    """
    Compute corse grained vorticity
    """
    corse = np.arange(50)
    
    for c in corse:
        data,field = vgradient.measure(M,'omega',step=c)
        name = 'omega_c'+str(int(c))
        print(name)
        setattr(M,name,data)
    
    
def movie(M):
    graphes.movie(M,'omega',Dirname=savedir)
    graphes.movie(M,'E',Dirname=savedir)
        
def track(M,sigma=3.):    
    pos = track.position(M,sigma=sigma)
    pos['sigma'] = sigma
    setattr(M,'core_pos',pos)

def iterate(Mlist,fun,*args,**kwargs):
    for M in Mlist:
        if M is not None:
            fun(M,*args,**kwargs)
        
def main(date,indices):
    Mlist = load(date,indices)    

    iterate(Mlist,compute)
   
    movie(Mlist[0])
    for M in Mlist:
        print(M)
        if M is not None:
            M.write(data=True,overwrite=True)
    
    

main(date,indices)