# -*- coding: utf-8 -*-

import glob
import os.path

"""
Module for exploring folders, save files and evaluate string (for title files typically)
"""

def get_fileList(fileDir,frmt,root='',display=False,sort='date'):
    fileName = fileDir + root+'*.' + frmt
    fileList = glob.glob(fileName)
    
    if sort=='name':
        fileList=sortbyname(fileList)
    else:
        fileList=sortbydate(fileList)
    
    n=len(fileList)
    if display:
        print('Number of ' + frmt + ' files : ' + str(n))
        print(fileName)
        
    return fileList,n

#get a number between two given separators : the first one should be unique in the string
def get_number(s,start,end,display=False,shift=0,from_end=False):
    val = get_string(s,start,end=end,display=display,shift=shift,from_end=from_end)
    if val=='':
        return -1
        
    try:
        num = float(val)
        return num
    except:
        print("cannot convert to float, keep it in the original format")
        return val
        
def get_end(s,start,end=0,shift=0,display=False):
    n=len(start)
    return s[n:]
    
def get_string(s,start,end='',shift=0,display=False,from_end=False):
    if from_end:
        j_p=str.find(s,end)
        i_p=str.find(s[j_p-1:0:-1],start[::-1])
                
        if j_p==-1:
            i=-1
            j=-1
            return ''
        else:
#            i=j_p-i_p-1
#            j = j_p-i-len(start)           
            spart = s[j_p-1:0:-1]
            return spart[i_p-1::-1]
    #    print(s[i+len(start):i+j+1])
    else:
        i=str.find(s,start)
      #  print(i)
      #  print(s)
        j=str.find(s[i+len(start):],end)
    
    if (i==-1)or(j==-1):
    #    print('failed')
        return ''
    else:
        deb=i+len(start)+shift
        fin=i+j+len(start)

        if display:
         #   print(deb)
         #   print(fin)
            print(s[deb:fin])

        return s[deb:fin]

def contain(s,pattern):
    return str.find(s,pattern)>=0
        
def get_dataDir(cinefile,root='/PIV_',frmat='cine',display=True):
    """
    Localize the data files associated to a cine file
        Look inside the directory and its sub-directory to find standard name of Data file folder.
        Several notation convention has been used, so this program look at all the different possibilities
    
    INPUT 
        cinefile : string, filename 
            name of the cine file
        root : depreciated (optionnal)
        frmat :  format of the initial data file (optionnal) default value is cine
        
    OUTPUT
        dataDirList : list of string. 
            A List of the folder names containing the PIV data
        n : number of folders
    """
    #loca
    #rootname of the cinefile : (without the extension and the PIV_ beginning)
    baseDir = os.path.dirname(cinefile)
    basename=get_string(os.path.basename(cinefile),'','.'+frmat)

    folders=['/PIV_data','/PIV_data_full','/PIV_step1_data','/PIV_step10_data','/PIV_data_extract']

    #root for PIV lab data. the folder have to start by the explicit keyword "PIVlab"
  #  print(basename)
    root='PIVlab'+'*'+basename
    
    dataDirList = []
    n=0
    
    for folder in folders:
        print("")
        print("")
        print("")
        print(folder)
        dataDirL = []
        Dir = baseDir + folder
        print('Looking for data in : '+Dir)
        print("toto")
        if os.path.isdir(Dir):
            if display:
                print("PIV data found at : "+Dir)
            dataDirL,n=get_dirList(Dir,root,display)

            #check if the associated dir is not empty
            keep_indices=[]
            for i,dataDir in enumerate(dataDirL):
                print("test : "+dataDir)
                
                f,n=get_fileList(dataDir+'/',frmt='txt')
             #   print(str(f))
                if not n==0:
                    print('ok !')
                    keep_indices.append(i)
            dataDirL = [dataDirL[i] for i in keep_indices]
            
        dataDirList += dataDirL
        
    n = len(dataDirList)
    if True:
        print("Datadir found :"+str(dataDirList))
    #find the list of file in the folder or the subfolder that match this pattern
    #look either in the same directory, or in a subdirectory named from the cinefile name
#    Dir=get_string(cinefile,'',root)    
#    print("Look for data inside : "+Dir)
#    if os.path.isdir(Dir):
#        print('looking directly in the cine directory')
#        Dir=os.path.dirname(cinefile)
#        root='PIVlab'+'*'+cineroot
#        dataDirList,n=get_dirList(Dir,root)
#    else:
#        print('looking inside the subdirectory')
 #       print(Dir)
#        dataDirList,n=get_dirList(Dir,'data_'+cineroot)
       # dataDirList,n=get_dirList(Dir,cineroot)
    return dataDirList,n

def get_dirList(Dir,root,display=False):
    #return a list of Directory contained in Dir, and sharing the same string root
    #dirList is by default sorted by date of creation (easier to index)
    fileParent=Dir + '/*'+root
    
    print(fileParent)
    DirParent=glob.glob(fileParent)
    
    
    fileSons=Dir + '/*'+root +'_*'
    DirSons=glob.glob(fileSons)
    DirList=DirParent+DirSons
    
    #How to sort them ?? -> also with the date of creation (unique)
    DirList=sortbydate(DirList)
    
    n=len(DirList)
    if display:
        print('Number of ' + root + ' Dir : ' + str(n))
        
    return DirList,n

def sortbydate(fileList):
    #current way of naming data files data_ + root_cinefile + others parameters if so
    #with a function calling directly the PIV measurements

    #sort fileList by date of creation
    tuple_dated=[(name, os.path.getctime(name)) for name in fileList]
    #sort by date of creation
    tuple_dated=sorted(tuple_dated,key=lambda t:t[1])
    fileList=[t[0] for t in tuple_dated]
    
    return fileList
    
def sortbyname(fileList):
    fileList.sort()
    return fileList
    
def digit_to_s(number,ndigit):
    s=str(number)
    
    for i in range(ndigit-len(s)):
        s='0'+s
    
#    print(s)
    return s
    
    #a database must be generate each day, to inventory all the previous experiments indexed in such a way