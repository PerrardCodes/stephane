import numpy as np
import matplotlib.pyplot as plt
import glob
import math

import stephane.jhtd.strain_tensor as strain
import stephane.jhtd.get as jhtd_get
import stephane.tools.dict2list as dict2list
import stephane.display.graphes as graphes

def process(data):
    """
    
    """
    Etot = jhtd_get.vlist(data)
    eigen,omega,cosine  = strain.strain_distribution(data)
        #data correspond to a slice in time, with nx x ny x nz points        
	#Total energy
  #      E=np.sum(np.power(data,2),axis=3)        
#    print('done')
    return Etot,eigen,omega,cosine

def energy_distribution(data):
    #E+=np.ndarray.tolist(np.power(data[key],2,axis=3))
    #plt.savefig(filename,dpi=300)
    filename = 'distribution_E_1.png'
    return []
        

def hist(Y,label,step=100,fignum=1):
    """
    Compute and plot the 1d histogramm in semilogy
    INPUT
    -----	
    Y : input data, list format
        in any format that can be translated by numpy.asarray
    label : string
        label to use for the plot
    step : relative number of histogramm steps (optional)
        step = 100 corresponds to N/100 bins where N is the total amount of data    
    OUTPUT
    ------
    None
    """    
    N=len(Y)
    n,x=np.histogram(np.asarray(Y),int(N/step))
    
    #normalize n :
    dx = np.mean(np.diff(x))
    n=n/(np.sum(n)*dx)
    plt.figure(fignum)
    
    plt.plot(x[1:],n,label)
    plt.axis([-1,1,0,np.max(n)*1.1])
    plt.pause(0.001)
   # plt.draw()
    
    plt.hold(True)
    plt.show(block=False)
        
def subset(eigen,omega,cosine,condition,val,type='>'):
    """
    Select a subset of the data corresponding to a given condition
        condition corresponds to the key of the conditionning variable
    INPUT
    -----
    eigen :
    omega :
    cosine :
    condition : float
    type : unused now
    OUTPUT
    -----
    """
    #print(cosine_t['lambda_omega_1'])
    indices = np.where(np.abs(np.asarray(cosine[condition]))>val)[0]
    print('Subset ratio : '+str(len(indices)*100./len(cosine['lambda_omega_2'])))

    eigen_extract = dict2list.extract_dict(eigen,indices)
    cosine_extract = dict2list.extract_dict(cosine,indices)
    
    #omega is given in another format, subset indices does not match for now
    #TO BE FIXED
    omega_extract = np.asarray(omega)#[indices]
    
    return eigen_extract,omega_extract,cosine_extract


 
def average_values(eigen,omega,cosine,d=3):
    """
    Compute the average values of strain tensor diagonal compnonents and enstrophy production 
    Can be compare qualitatively to Tsinober 1998 (both experiments and old DNS data)
    outputs are printed in the console
    INPUT
    -----
    eigen : Dictionnary containing the eigenvalues Lambda and eigenvectors lambda
    omega : Dictionnary containing components of the vorticity field
    cosine : Dictionnary containing orientation angle between lambda and omega
    OUTPUT
    -----
    None
    """
    omega_2 = np.asarray(dict2list.to_1d_list(np.power(omega,2,),d=4))
#    omega_2 = np.asarray(to_1d_list(np.sum(np.power(omega,2),axis=4),d=4))

    #Even quantities
    values = [0.32, 0.15, 0.53]
    val=np.zeros(3)
    
    print(eigen.keys())
    for i in range(d):
        var=np.power(eigen['Lambda_'+str(i)],2)*omega_2[:,i]*np.power(cosine['lambda_omega_'+str(i)],2)
  #      print(var.shape)
        val[i] = np.mean(var)

    val = val / np.sum(np.power(val,1))
    
    print("")
    print("<omega^2 L^2 cos^2 > :")
    for i in range(d):
        print(str(val[i]) +" ("+str(values[i])+")")
    
    values = [0.38, 0.11, 0.51]    
    print("")
    print("<L^2cos^2 > :")
    for i in range(d):
        Lambda_i = eigen['Lambda_'+str(i)]
        cos_i = cosine['lambda_omega_'+str(i)]
        val[i] = np.mean(np.power(Lambda_i,2)*np.power(cos_i,2))
        
    val = val / np.sum(np.power(val,1))
    for i in range(d):
        print(str(val[i]) +" ("+str(values[i])+")")
    
    #odd quantities
    values = [-0.63, 0.46, 1.17]    
    print("")
    print("<L_i cos^2 > :")
    for i in range(d):
        Lambda_i = eigen['Lambda_'+str(i)]
        cos_i = cosine['lambda_omega_'+str(i)]
        val[i] = np.mean(Lambda_i*np.power(cos_i,2))
        
    val = val / np.sum(np.power(val,1))
    for i in range(d):
        print(str(val[i]) +" ("+str(values[i])+")")
    
    
    values = [-0.51, 0.51, 1.06]    
    print("")
    print("<omega^2 L_i cos^2 > :")
    for i in range(d):
        Lambda_i = eigen['Lambda_'+str(i)]
        cos_i = cosine['lambda_omega_'+str(i)]
        val[i] = np.mean(Lambda_i*omega_2[:,i]*np.power(cos_i,2))
        
    val = val / np.sum(np.power(val,1))
    for i in range(d):
        print(str(val[i]) +" ("+str(values[i])+")")

def norm(u,axis=1):
    return np.sqrt(np.sum(np.asarray(u)**2,axis=axis))


def plots(eigen,omega,cosine,step):
    """
    Make plots of geometrical quantities associated to the strain tensor 
    (eigenvalues, vorticity and stretching vector)
    INPUT
    -----
    eigen : Dictionnary containing the eigenvalues Lambda and eigenvectors lambda
    omega : Dictionnary containing components of the vorticity field
    cosine : orientation angle between lambda and omega
    step : average number of data point per bin
    OUTPUT
    -----
    figs : dict
        dictionnary of output figures, the key correspond to the number of the figure
        associated value is a title in string format (root name for an eventual saving process)
    """
    figs={}
    
    #print('Epsilon : ')
    graphes.hist(eigen['epsilon'],label='k',step=step,fignum=1)
    figs.update(graphes.legende('$\epsilon$','PDF','',display=False))
    
    label = ['k','b','r']
    if True:
#    for i,key in enumerate(eigen.keys()): 
      #  k='Lambda_'
      #  if key.find(k)>=0:
      #      j=int(key[len(k)])
            #hist(eigen_t[key],label=label[j],step=step,fignum=2+j)
            #plt.title(key)
        enstrophy = norm(omega,axis=3)
        graphes.hist(enstrophy,label='r',step=step,fignum=2)
        figs.update(graphes.legende('$\omega$','PDF','',display=False))
        
#if False:    
        for i,key in enumerate(cosine.keys()): 
        #    print(key)
            keys = ['lambda_omega_','lambda_W_']
    
            for z,k in enumerate(keys):
                if key.find(k)>=0:
                    j=int(key[len(k)])
                 #   print(j)
                    graphes.hist(cosine[key],label=label[j],step=step,fignum=5+3*z+j)
                    if z==0:
                        figs.update(graphes.legende('cos($\lambda_'+str(3-j)+',\omega$)','PDF','',display=False))
                    if z==1:
                        figs.update(graphes.legende('cos($\lambda_'+str(3-j)+',W$)','PDF','',display=False))    
                                    
            if key.find('W_omega')>=0:
 #               print(step)
                graphes.hist(cosine[key],label='k',step=step,fignum=15)   
                figs.update(graphes.legende('cos($\omega,W$)','PDF','',display=False))
    
#    print(figs)       
    return figs

if __name__== '__main__':
        """
        load the data file by file, and compute vorticity, eigenvectors of the strain tensor and geometrical angles for each single point
        Call plots() and average_values() as outputs
        """
    #for i in range(1):
    #    file ='/Users/stephane/Documents/JHT_Database/Data_sample_2015_10_07/Tests_isotropic1024coarse_100_15_15_15_'+str(i)+'.h5'
        Datadir = '/Volumes/labshared/JHTDB_data/JHTDB_lib/'#/Users/stephane/Documents/JHT_Database/Programs_JH/JHTDB_lib'+'/'
    
        fileList = glob.glob(Datadir+'*.h5')
        N = len(fileList)

        E_t=[]
        omega_t=[]

        cosine_t={}
        eigen_t={}

        count=0

        N=24000
        step = 500

        print(len(fileList))
        for i,file in enumerate(fileList[:N]):

            if i%100==0:
                print(str(i*100//N)+' %')
    #    file ='/Users/stephane/Documents/JHT_Database/Data_sample_2015_10_07/Serie/isotropic1024coarse_'+str(i)+'.h5'  
            try:
                data = read_jhdt(file)   
            except:
                count+=1
           # print("Datafile not opened : "+file)
    
            E,eigen,omega,cosine = process(data)

            E_t += E    
            omega_t += omega 
    
            eigen_t = dict2list.add_list_to_dict(eigen_t,eigen)
            cosine_t = dict2list.add_list_to_dict(cosine_t,cosine)  
    

        print("Percentage of data files corrupted : "+str((100*count)/len(fileList))+" %")

        plots(eigen,omega,cosine,step)

        T=0.9
    #print(cosine_t['lambda_omega_1'])
        indices = np.where(np.abs(np.asarray(cosine_t['lambda_omega_1']))>T)[0]

        print(len(indices)*100./len(cosine_t['lambda_omega_2']))

        eigen_extract = dict2list.extract_dict(eigen_t,indices)
        cosine_extract = dict2list.extract_dict(cosine_t,indices)
        omega_extract = omega_t#np.asarray(omega_t)[indices]

        plots(eigen_extract,omega_extract,cosine_extract,step)

        average_values(eigen,omega,cosine)