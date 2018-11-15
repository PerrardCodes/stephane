


import numpy as np
import stephane.display.graphes as graphes
import matplotlib.pyplot as plt
import stephane.tools.dict2list as dict2list


def profile_1d(M,Dt=10,direction='v',start=20,fignum=1):
    """
    Compute the 1d averaged profile> Averaging is performed along a Dt window in time, and along the specified axis, either 'v' or 'h'
    INPUT
    -----
    M : Mdata object.
    Dt : int. Default 10
        time windows width for the averaging
    direction : string. default 'v'. Only option now
    start : int
        starting index
    fignum : int
        number for the figure output
    OUTPUT
    -----
    Ux,Uy,Ux_std,Uy_std
    """
    dimensions=M.shape()

    if M.param.typeplane=='vp':
        if M.param.angle==0.:
            z=M.y[:,0]
            axis = [1,0]
        if M.param.angle==90.:
            z=M.x[0,:]
            axis = [0,1]
    else:
        z=M.y[:,0]
        axis = [1,0]
            
    for i in range(start,dimensions[2],Dt):
        #averaging over Dt in time, and along one dimension in space
        Ux=np.nanmean(np.nanmean(M.Ux[...,i:i+Dt],axis=2),axis=axis[0])
        Uy=np.nanmean(np.nanmean(M.Uy[...,i:i+Dt],axis=2),axis=axis[0])
        
        #standard deviation computation
       # print(dimensions[axis[1]])
    #    print(tuple(axis+[Dt]))
        
        orientation = tuple([axis[0]+1,axis[1]+1]+[0])
        Ux_mean = np.asarray(np.transpose([[dict2list.to_1d_list(Ux) for k in range(dimensions[axis[0]])] for t in range(Dt)],orientation))
        Uy_mean = np.asarray(np.transpose([[dict2list.to_1d_list(Uy) for k in range(dimensions[axis[0]])] for t in range(Dt)],orientation))

#        Uy_mean = np.asarray([[dict2list.to_1d_list(Uy) for k in range(dimensions[axis[0]])] for t in range(Dt)],tuple(axis+[Dt]))        
        std_Ux=np.sqrt(np.mean(np.mean(np.abs(M.Ux[...,i:i+Dt]-Ux_mean)**2,axis=2),axis=axis[0]))
        std_Uy=np.sqrt(np.mean(np.mean(np.abs(M.Uy[...,i:i+Dt]-Uy_mean)**2,axis=2),axis=axis[0]))
  
        graphes.set_fig(0) #clear current axis
        graphes.graph(z,Ux,label='k^',fignum=fignum)#,std_Ux)
        graphes.graph(z,Uy,label='ro',fignum=fignum)#,std_Ux)
        graphes.set_axis(-400,-100,-100,100)
        
        fig = graphes.legende('z (mm)','V_{rms} (mm/s)','Ux, Uy')
        
        filename = './Results/Mean_profile_'+M.id.get_id()+'/'+ fig[fignum] + '_t_' + str(i)
       # print(filename)
        
        graphes.save_fig(fignum,filename,frmt='png')
#        graphes.graph(z,Uy,std_Uy)
#        graphes.legende('z (m)','V (m/s)','Uy')
      #  raw_input()
        
        
    return Ux,Uy,std_Ux,std_Uy
    
    
def time_average(M,i,Dt=50):
    """
    Compute the profile averaged over time
    """
    

def ensemble_average(Mlist,i,Dt=10):
    """
    Compute an ensemble average from a set of identical experiments
        for now, superposed the fields (require exactly the same set-ups)
    """
    
    
    