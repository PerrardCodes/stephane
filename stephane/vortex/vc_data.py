import stephane.vortex as vortex


"""
Generate a vc_data type of data (stored in a dictionnary, and make a local folder storing all these trajectories)

"""


def generate(date,indices):
    
    Slist = Sdata_manip.load_serie(date,indices=[0])
    Mlist = Sdata_manip.load_measures(Slist,indices=0)
    
    for M in Mlist:
        M.get('E')
        M.get('omega')
    
    for M in Mlist:
        if hasattr(M,'Id'):
            if M.Id=='pointer':
                M.Id = M.Sdata.Id
    pos_total = []
    tmin = 20
    tmax = 158

    for i,M in enumerate(Mlist):
        pos = track.position(M,sigma=2.)
        pos['param'] = add_parameters(M)


        accurate=[]
        for i,M in enumerate(Mlist):
            figs,acc = track.plot(pos[i],tmin,tmax)
            accurate.append(acc)
        graphes.save_figs(figs,savedir=savedir+'')
        
        
def add_parameters(M):
    Dict = {}

    Dict['cinefile'] = M.Sdata.fileCine
    Dict['date'] = M.Id.date
    Dict['index'] = M.Id.index
    Dict['A'] = browse.get_number(Dict['cinefile'],'_A','mm')
    Dict['v'] = browse.get_number(Dict['cinefile'],'_v','mms')
    if Dict['v']==-1:
        v = browse.get_number(cinefile,'mm_v','mm')
    Dict['fps'] = browse.get_number(cinefile,'_fps','_')
    Dict['fx'] = 

    
    