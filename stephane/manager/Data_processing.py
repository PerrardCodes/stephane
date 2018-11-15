# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 11:20:23 2015
Generates :
the timestep files for sampling the images
generate the tiff files following the instruction given in timestep
@author: stephane
"""


#Process data
#several steps to process data :
# generates ref movies for each cine file in the date folder

import stephane.manager.ref_movie as ref_movie
import stephane.manager.file_architecture as file_architecture
import stephane.manager.cine2pic as cine2pic

import stephane.tools.browse as browse
import stephane.tools.rw_data as rw_data
import stephane.cine as cine

import argparse
import os.path

import glob

#script_dir = os.path.dirname(os.path.realpath(__file__))
parser = argparse.ArgumentParser(description="Process PIV data date by date")
parser.add_argument('-d',dest='date',default='1987_03_19',type=str,help='date to be processed. Python will look for it in the folders specified in file_architecture.py')
parser.add_argument('-f',dest='folder',default=None,type=str,help='base folder to be processed. Python will look for cine files inside this specified folder')
parser.add_argument('-s',dest='start',default=None,type=int, help='start processing index of the cinefile List')
parser.add_argument('-e',dest='end',default=None,type=int, help='end processing index of the cinefile List')
parser.add_argument('-n',dest='n',default=10,type=int, help='Number of images for the ref movie (optionnal)')
parser.add_argument('-step',dest='step',default=2,type=int, help='Under sampling of the data. Default value is 2')
args = parser.parse_args()

def main(date):
    """
    main function
    INTPUT
    -----
    date : string. 
        Date of the experiments to be processed. 
        The associated folder is researched by file_architecture.py where the possible locations of the data on immense are listed
    OUTPUT
    -----
    None    
    """
#    date = args.date
    #if the directory is not specified, look automatticlly in grid_turbulence folders
    if args.folder is None:
        args.folder = file_architecture.get_dir(date)
    process(args.folder)
    list_item = glob.glob(args.folder+'/*')
    
    #Recursive search
    for l in list_item:
        if os.path.isdir(l):
            process(l)
        
def process(directory):
    """
    Generate tiff image folders for cine files located in a specified directory
    INPUT
    -----
    folder : string
        path of the folder to process
    OUTPUT
    -----
    """
    fileList,n = browse.get_fileList(directory,frmt='cine',root='/PIV_')
    fileList_ref,n = browse.get_fileList(directory,frmt='cine',root='/Reference_')
    fileList_bubble,n = browse.get_fileList(directory,frmt='cine',root='/Bubble')
    fileList_other,n = browse.get_fileList(directory,frmt='cine',root='/')
    
    print(fileList_other)
  #  print(fileList)
    #print(fileList_ref)
    
    step=0
    while step>=0:
        if step == 0:
            ref_movie.make_ref(fileList_ref)
            ref_movie.make_ref(fileList_bubble)
#            caller(fileList_ref,step)
    #    step = caller(fileList,step)
        step = caller(fileList_other,step)
        

        
def caller(fileList,step):
    """
    call the nth function of a function list
    INPUT
    -----
    fileList : list of cine file name to be processed
    step : int. Number of the function to call
    -----
    OUTPUT
    step : updated value of step. ++1 if still running, -1 if all the functions have been called
    """
    list_fun = [ref_movie.make,         # generate a movie directory to be used as a reference (grid position, spatial scale, eventually time zero)
                make_timestep_files,    # generate timestep txt files to be used as an input for generating to-be-PIV-processed images
                make_piv_folders]#,       # input folder 
                #make_result_folder]     # output folders for PIV measurements. Everything is stored into a PIV_data file
                                       
    n = len(list_fun)
    functions = {key:fun for key,fun in zip(range(n),list_fun)}
                     
    if step<n:
        print(functions[step].__name__)
        functions[step](fileList)
        step+=1
    else:
        step=-1
    return step
    
def iterate(fileList,function):
    for file in fileList[args.start:args.end]:
        function(file)

def make_timestep_files(fileList,Dt=1,starts=[],ends=[],Dt_list=[]):
    """
    Generate a timestep .txt file that will be used to process the data with the right time step
    INPUT
    -----
    fileList : list of filename to process
    Dt : int. timestep to be applied. Default value : 1
        Rmq : could be switched to a list of timestep and associated start/end indexes for each timestep 
    Optionnal variables (not implemented yet) : starts, ends and Dt_list to set a list of timestep for different instants in the movie.
    OUTPUT
    -----
    NONE
    """
    
    keys = ['start','end','Dt']
    
    for file in fileList[args.start:args.end]:
        
        n = str(browse.get_string(file,start = '_n',end='_'))
        try:
            int(n)
        except:
            n=''
        if n=='':
            #get the number of images from the cinefile
            c=cine.Cine(file)
            n=c.len()-1
            print(n)
               
        values = [[0],[n],[Dt]]
        
        base = browse.get_string(os.path.basename(file),'','.cine')
        file_timestep = os.path.dirname(file) +'/PIV_timestep_'+base+'.txt'
        
        if not os.path.isfile(file_timestep):
            rw_data.write_dictionnary(file_timestep,keys,values)
     
def make_piv_folders(fileList,step=None):
    """
    Generate tiff folders associated to each cine file. 
        the timestep is given by the associated timestep.txt file previously generated
    INPUT
    -----
    fileList :
    step : int.  Default value is 20
        Decimation parameter, in order to save only 2/step of the images.
        (default value is 2 as two images are saved for each number (images A and B), spaced by the specified Dt)
    OUTPUT
    -----
    None
    """
    
    if step==None:
        step=args.step
    #from a dirname corresponding to a particular date, 
    #generate folders containing individual TIFF for PIV processing (no matter the PIV software used)
    for file in fileList[args.start:args.end]:
        print(file)
        cine2pic.cine2tiff(file,'File',step,post='_File')    
    print(str(len(fileList))+' cine files to be processed')

def make_result_folder(fileList,type_analysis='Sample',algo='PIVlab',W=32,ratio=None):
    """
    Create folders for storing the PIVlab Data (now processed in matlab)    
    INPUT
    -----	
    fileList : List of string
        filename List to be read
    type_analysis : string (default value 'Sample')
        Standard presaved types of analysis (Sample = 1/10 images processed, Full_single', every images, Full_double every pair of image (effective ratio = 2))
                Full_single : ratio = 1
                Full_double : ratio = 2
                Sample : ratio = 10
    algo : string (default value 'PIVlab')
        Name of the PIV algorithm used to processed the data
    W : int (default value = 16)
        Window size used for PIV processing.
    ratio : number of images / number of pair of images processed. default value = 10. Can be set to standard values using type_analysis
    OUTPUT
    ------
    None
    0	"""    
    types = dict(Sample=ratio,Full_double=2,Full_single=1)
    if type_analysis in types.keys():
        ratio = types[type_analysis]
    else:
        #use the value given in argument for ratio
        pass
    
    for file in fileList[args.start:args.end]:
        Dir = os.path.dirname(file)
        rootDir = Dir +'/PIV_data'
        basename = browse.get_string(os.path.basename(file),'',end='.cine')
        foldername = rootDir+'/'+algo+'_ratio'+str(ratio)+'_W'+str(W)+'pix_'+basename

        if not os.path.isdir(foldername):
            os.makedirs(foldername)
        
    #default    
    
#date = '2015_09_28'
#main(date)    
if __name__ == '__main__':
    main(args.date)
# By hand, measure the parameters associated to the cine file. Generate a .txt parameter file for each movie
# The ref files could bve generated automatically with 0 values (?!)

# Write Sdata associated to each of its cine file

# Launch PIV on each Sdata ? -> for now this step is computed outside of the pipeline in matlab
# Generate Mdata from the different PIV folders