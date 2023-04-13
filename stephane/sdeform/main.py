from stephane.sdeform.spherical import *

import numpy as np
import glob

from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
import matplotlib.pyplot as plt
import numpy as np
import stephane.tools.rw_data as rw_data



def main(folder=None):
    if folder == None:
        folder = '/Users/stephane/Dropbox/ARiviere/graphes/Bulle_vs_time/'

    filelist = glob.glob(folder +'*.dat')
    print(filelist)
    
    
    D = read_bubbledata(filelist)

    plt.hist(D['volume'])
    V0 = 1000
    plt.vlines(1000,ymin=0,ymax=16000)
    plt.show()
    
    
    # remove small bubbles from the data by giving a threshold in volume
    indices = np.where(np.asarray(D['volume'])>V0)
    D['indices'] = indices
    coords = ['x','y','z']
    for key in coords:
        D[key+'c'] = np.asarray(D[key+'0'])[indices]
    D['V']=np.asarray(D['volume'])[indices]
    
