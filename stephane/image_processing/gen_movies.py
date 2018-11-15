


import subprocess as sub
import glob
import os
import os.path
import stephane.Image_processing.cine2avi as cine2avi

def compress_movie(filename,rate = 4096):
    """
    Compress movie using ffmpeg. About a factor 30 of compression for default rate. Can be scale up or down on demand
    INPUT
    -----
    filename : str
        full filename of the avi file.
    rate : int 
        rate of data in kbit/s. default value is 4096, corresponding to a compression of about a factor 30.
    MORE OPTIONS TO COME !
    Warning : in case of overwriting, it will pop up a message in the terminal which stop the processing
    OUTPUT
    -----
    None
    
    """
    function = ['ffmpeg','-i']
    parameters = ['-c:v','libx264','-pix_fmt','yuv420p','-profile:v','high', '-level:v','5.0','-b:v',str(rate)+'k']

    output_name,ext = os.path.splitext(filename)
    output_name = output_name + '.mp4'
    
    sub.call(function+[filename]+parameters+[output_name])
        
def folder(dirname):
    """
    Generate a mp4 movie for each avi file contained in the given directory
    INPUT
    -----
    dirname : str
        path of a directory
    OUTPUT
    -----
    None
    """
    fileList = glob.glob(dirname+'/*.avi')
    print('Avi files found : '+str(fileList))
    for filename in fileList:
        compress_movie(filename)

def rm_avi(dirname):
    """
    Delete all the avi files contained in the given directory
    INPUT
    -----
    dirname : str
        path of a directory
    OUTPUT
    -----
    None
    """
    fileList = glob.glob(dirname+'/*.avi')
    print('Remove avi files :')
    print(fileList)
    
    for filename in fileList:
        os.remove(filename)
    
    print('done')
#def main():
#    folder(dirname)

def make(dirname,type='multiple',rect='[150:650:0:1280]',timestamp=True,quality=100,**kwargs):
    """
    Generate avi and mp4 movies for all the cine files contained in the given directory.
    Can be either in multiple modes (stack of movies for each final avi file) or single mode (one avi per each cine file) 
    INPUT
    -----
    dirname : str
        path of a directory
    type : str
        either 'single' or 'multiple'. default 'multiple'
    rect : str
        rectangle to crop the image. Command line format ([x0:x1:y0:y1]). default (current) '[150:650:0:1280]'
    timestamp : bool
        add a timestamp on top of the movie. default True
    quality :  int
        quality from 0 to 100. default 100
    **kwargs : not yet implemented
        add other parameters for generating avi (see cine2avi for details)
    OUTPUT
    -----
    None
    """
    fileList = glob.glob(dirname+'/*.cine')
    kwargs = dict(type=type,rect=rect,timestamp=timestamp,quality=quality)
    
    if type=='multiple':
        #only make avi files for series
        fileList = glob.glob(dirname+'/*_1.cine')
        #processed = []
        for filename in fileList:    
            base,ext = os.path.splitext(filename)
            base = base[:-1]+'*'
            filename = base + ext
            #if not filename in processed:
            print(filename)
            #processed.append(filename)
            try:
                cine2avi.make(filename,**kwargs)    
            except:
                print('Error in the avi generation file')
                
            #convert into mp4
            folder(dirname)
            #delete original avi movies
            #rm_avi(dirname)
    else:
        for filename in fileList:
            try:
                cine2avi.make(filename,**kwargs)    
            except:
                print('Error in the avi generation file')
            #convert into mp4
            folder(dirname)
            #delete original avi movies
            #rm_avi(dirname)  
              
def example_1(rootdir):
    """
    generate an avi file for each serie of cines contained in any subfolder of rootdir
    """
    dirnames = glob.glob(rootdir+'*')
    for dirname in dirnames:
        make(dirname,type='multiple',rect='[150:650:0:1280]',timestamp=True,quality=100,**kwargs)
    
def example_2(dirname):
    """
    generate an avi file for each cine. No cropping
    """
    make(dirname,type='single',rect=None,timestamp=True,quality=100)
    
def main():
    """"
    Currently designed to be launched directly from Irvinator-2.
    gen_movies can be imported in python and therfore use as a module. example_1() or example_2() can also be used
    """
    
    root = '/home/steph/Documents/Stephane/Vortex_Turbulence/'
    dirnames = glob.glob(root+'*')
    
    for dirname in dirnames:
        make(dirname)
         
if __name__ == '__main__':
    main()    