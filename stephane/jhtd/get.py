
import h5py
import numpy as np
import glob

import stephane.tools.rw_data as rw_data
import stephane.tools.browse as browse
import stephane.jhtd.cutout as cutout


def read(file,key_type = u'u'):
    """
    Read a JHDT file (h5py format) 
        read a h5py file
        extract from that the dataset names u'uX' where X denotes the time of the frame
    INPUT
    -----	
    file : string
        filename of the h5py file
    OUTPUT
    ------
    return a dictionnary with same keys as the initial set, each field containing one time step
    	"""
        
    f=h5py.File(file,'r')
   # print(f.keys())
    data_keys=[key for key in f.keys() if key.find(key_type)>=0]
    data={}
    for key in data_keys:
        data[key]=f[key]
    return data
    
def read_chunk(folder,date):
    """
    Read all the data contained in the subfolder of the called folder
    """
#    folder = '/Users/stephane/Documents/JHT_Database/Data/Spatial_measurement_1d_2016_03_30/Data'
    l = glob.glob(folder+'/zl*')
    print(folder)
    data = {}
    param = {}
    files_fail = []
    for i in range(len(l)):
       # print(i)
        files = glob.glob(l[i]+'/*.hdf5')
        for file in files:
            try:
                data_part = read(file)
                param[data_part.keys()[0]]=get_parameters(file)
                data.update(data_part)
                
            except:
                print(str(i) + ', skipped')
             #   print(file)
                files_fail.append(file)

    print('Number of failed loading : '+str(len(files_fail)))
   # print(files_fail)
    #log_file = generate_log(files_fail,date)    
    #print(log_file)
    
#    Header,dset_list = rw_data.read_dataFile(log_file,Hdelimiter='\t',Ddelimiter='\t',Oneline=False)

   # cutout.recover(log_file)
    return data,param

def generate_log(files_fail,date):
    dirbase = '/Users/stephane/Documents/JHT_Database/Data/Spatial_measurement_2d_'+date+'/'
    List_key = ['zl', 'yl', 'xl','t0','tl','y0','x0','z0']

    log_file = dirbase + 'log.txt';
    f_log = open(log_file,'w')
    rw_data.write_header(f_log,List_key) #write the header

    for file in files_fail:
    #    print(file)
        param = get_parameters(file)
        rw_data.write_line(f_log,param)
    f_log.close()
    print(log_file + ' generated')
    
    return log_file

def get_parameters(file):  
    List_key = ['zl', 'yl', 'xl','t0','tl','y0','x0','z0']
     
    param = {}
    for k in List_key[:-1]:
        param[k] = int(browse.get_string(file,'_'+k+'_',end='_',display=False))
    k=List_key[-1]
    param[k] = int(browse.get_string(file,'_'+k+'_',end='/',display=False))
    
    return param    
    
def vlist(data,rm_nan=True):
    """
    Extract the velocity components from a JHTD data format. Spatial coherence is lost (return a list of single data points)
    INPUT
    -----
        data : JHTD data format
    OUTPUT
    -----
        U : numpy array of 1d1c velocity components
             each element corresponds to one point 3C of velocity
    """
    U=[]
    for key in data.keys():
	    #wrap the three components in one 3d matrix (by arbitrary multiplying the first dimension by 3 : spatial                information is lost)
        dim=data[key].shape
        dim=(dim[0]*3,dim[1],dim[2],1)
        Upart=np.reshape(data[key],dim)

        #wrap to 1d vector every single component
        dimensions = Upart.shape
        N = np.product(np.asarray(dimensions))
        U_1d = np.reshape(Upart,(N,))
        if rm_nan:
            if np.isnan(U_1d).any():
                print('Nan value encountered : removed from the data')
            U_1d = U_1d[~np.isnan(U_1d)]
            
        U += np.ndarray.tolist(U_1d)

    return U
    
def generate_struct(data,param):
    
    JHTD_data = JHTDdata(data,param)
    return JHTD_data
    #from a dataset, generate a Mdata structure
    
class JHTDdata:
    def __init__(self,data,param,N=None):
        import stephane.mdata.Id as Id
        
        self.load_data(data,param,N=N)
     #   self.filename = 
        self.Id=Id.Id(S=None,typ='Numerics',who='JHTDB')

    def load_data(self,data,param,N=None):
        if N==None:
            N = len(data.keys())
        
        dim = data[data.keys()[0]].shape
        tup = tuple(dim[:-2])+(N,)
        
     #   print(tup)
        
        self.t = np.zeros(N)
        self.x0 = np.zeros(N)
        self.y0 = np.zeros(N)
        self.z0 = np.zeros(N)
        
        self.Ux = np.zeros(tup)
        self.Uy = np.zeros(tup)
        self.Uz = np.zeros(tup)

        self.def_axis(data,param,dim)
        
        self.fx=1
        self.ft=1
        
        keys = np.sort(data.keys())
        for i,key in enumerate(keys[:N]):
            self.t[i]=param[key]['t0']
            self.x0[i]=param[key]['x0']
            self.y0[i]=param[key]['y0']
            self.z0[i]=param[key]['z0']
            
            
            data_np = np.asarray(data[key])

            self.Ux[...,i] = data_np[...,0,2] 
            self.Uy[...,i] = data_np[...,0,1]
            self.Uz[...,i] = data_np[...,0,0]
        
        indices = np.argsort(self.t)
        
        self.t = self.t[indices]
        self.x0 = self.x0[indices]
        self.y0 = self.y0[indices]
        self.z0 = self.x0[indices]

        self.Ux = self.Ux[...,indices]
        self.Uy = self.Uy[...,indices]
        self.Uz = self.Uz[...,indices]        
        
    def def_axis(self,data,param,dim,d=2):
        key = data.keys()[0]
        #first dimension is x instead of z
        self.x=np.asarray([[[k for i in np.arange(param[key]['xl'])] for j in range(param[key]['yl'])] for k in range(param[key]['zl'])])[...,0]
     #   print(self.x)
        self.y=np.asarray([[[j for i in np.arange(param[key]['xl'])] for j in range(param[key]['yl'])] for k in range(param[key]['zl'])])[...,0]
        self.z=np.asarray([[[i for i in np.arange(param[key]['xl'])] for j in range(param[key]['yl'])] for k in range(param[key]['zl'])])[...,0]
                
     #   print(self.x.shape)
     #   print(self.y.shape)
     #   print(self.z.shape)

#        self.x = np.ones(tuple(dim[:-2]))
#        self.y = np.ones(tuple(dim[:-2]))
#        self.z = np.ones(tuple(dim[:-2]))
        
    def shape(self):
        return self.Ux.shape
        
        
    
########################## Measurements #############################

    def get(self,field,**kwargs):
        
        if field == 'U':
            #return both component in a vectorial format
            Ux = self.get('Ux')
            Uy = self.get('Uy')
            data = np.transpose(np.asarray([Ux,Uy]),(1,2,3,0))
            return data
#        if (not hasattr(self,field)) or (compute):
            
#            vgradient.compute(M,field,Dt=Dt_filt,**kwargs)        
        if not hasattr(self,field):
            if 'Dt_filt' in kwargs and kwargs['Dt_filt']>1:
                print('Filtering of the data : irreversible')
            self.compute(field)
#            setattr(self,field,)
        if hasattr(self,field):
            return getattr(self,field)
        else:
            return None
    
    def get_cut(self):
        print("crazy")
        
    def compute(self,field,**kwargs):
        import stephane.analysis.vgradient as vgradient
        return vgradient.compute(self,field,**kwargs)

    def measure(self,name,function,force=False,*args,**kwargs):
        if (not hasattr(self,name)) or force:
            print("Compute "+name)
            val = function(self,*args,**kwargs)
            setattr(self,name,val)
        else:
            print("Already computed")    
    