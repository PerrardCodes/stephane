# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 14:15:20 2015

@author: stephane
"""

import os.path
import numpy as np
#to use to read data files in ASCII formats (or others ?)
#to use to write data files in ASCII formats (or others ?) with label on the top :
# in particular to create a catalog of the parameters of all the existing data (!)

def gen_name(param_list):
    """
    generate a standart formatted filename from a dictionnary of parameters given in input
	
    INPUT
    -----	
    param_list : dictionnary of string convertible object
        List of parameter used
 
    OUTPUT
    ------
    Return a string generated from the dictionnary of parameters
    	"""

    string=''
    for key in param_list.keys():
        string=string+'_'+key+str(param_list[key])+'unit'
    return string

def read_dataFile(file,Hdelimiter=',',Ddelimiter=',',Oneline=False):
    """
    Basic function for reading a .txt dataFile
	
    The data file format has to be as follow :
    Header of one or more lines containing the title of each dataset in the last Header line.
    Data are presented into columns, separeted by a given character Ddelimiter
    
    INPUT
    -----	
    file : string
        filename to be read
    Hdelimiter : string (default value ',')
        Delimiter used between two adjacent values inside the Header. Usual values ',' '\t'
    Ddelimiter : string (default value ',')
        Delimiter used between two adjacent values inside the Data. Usual values ',' '\t'
    Oneline : bool (default value False)
        If Oneline, then return single elements instead of lists (is it a good idea ??)

    OUTPUT
    ------
    Return a tuple containing the Header of the file, and Data in a dictionnary format. Each key is given by the title present at the top of each column
    	"""    
    f = open(file,'r')
    Header,fline=read_begin(f,Hdelimiter)
    #Name the fields from the last line of the Header
    Names=parseline(Header[-1],Hdelimiter)
   # print(Names)
    #Generate a dictionnary with these names
   # Names=[name[0] for name in Names]

    Names=[name for name in Names]    
    Data={name:[] for name in Names}  #Just keep the first letter for now (!) should be improved for the first list of characters without space
        
    Data=read_line(fline,Names,Data,Ddelimiter) 
    for line in f:
        Data=read_line(line,Names,Data,Ddelimiter)
        
    f.close()
    return Header,Data
    
def read_matrix(file,Hdelimiter=',',Ddelimiter=','):
    """
    Read a datafile containing two entries data
	
    The data file format has to be as follow :
    Header of one or more lines containing the title of each dataset in the last Header line.
    Data are presented into columns, separeted by a given character Ddelimiter.
    Two of the data line correspond to both axis of the matrix, so that they contain equal values (respectively nx and ny similar values)
    
    INPUT
    -----	
    file : string
        filename to be read
    Hdelimiter : string (default value ',')
        Delimiter used between two adjacent values inside the Header. Usual values ',' '\t'
    Ddelimiter : string (default value ',')
        Delimiter used between two adjacent values inside the Data. Usual values ',' '\t'

    OUTPUT
    ------
    Return a tuple containing the Header of the file, the Data in a dictionnary format, and the two detected axis that present nx x ny equal values.
    	"""    
    Header,Data = read_dataFile(file,Hdelimiter=Hdelimiter,Ddelimiter=Ddelimiter)
    
    #bad way of reading file : a header should be added on each text file with the information needed to unwrap the matrix
    #Unwrap the Data to express it as matrices.    
    axes={}
    for i,key in enumerate(Data.keys()):
        count = 0
        for d in Data[key]:
            if d == Data[key][0]:            
                count=count+1
        if count > 1:
            N = len(Data[key])
            axes[key] = N//count
            print(key +' : '+str(axes[key]))

    return Header,Data,axes
    
            
def read_Header(file,delimiter):
    """
    Read the Header of a file.
    
    The Header is defined recursively by non numerical values. Once a numerical value is found, the Header end is set to the previous line. 
    
    INPUT
    -----	
    file : string
        filename to be read
    delimiter : string (default value ',')
        Delimiter used between two adjacent values inside the Header. Usual values ',' '\t'

    OUTPUT
    ------
    Return a string of the entire Header
    	""" 
    f = open(file,'r')    
    Header,fline=read_begin(f,delimiter)
    f.close()
    return Header
    
def read_begin(f,delimiter):
    #load the first lines, while not finding number in the first line element (fline element should not be a number in the Header !)
    #when a number is finally found,m return the Header and the first line of values
    number=False
    Header=[]
    while not number:
        fline=f.readline()
        List=parseline(fline,delimiter)
        try:
            #try only for the first splitted element
           # print(List[0])
            float(List[0])
            number=True
        except ValueError:
           # print(fline)
            #print("Not a number, pass to the next line")
            number=False
            Header.append(fline)
    
    return Header,fline
    
def read_line(line,Names,Data,delimiter):
    #require float everywhere  !
    List=parseline(line,delimiter)
    
#    print(Data)
#    print(Names)
    for i,v in enumerate(Data):
        if Names.index(v)<len(List):
            try:
                Data[v].append(float(List[Names.index(v)]))
            except:
                #print('could not convert '+str(List[Names.index(v)])+' to float')
                Data[v].append(List[Names.index(v)])
        else:
        #    print(v)
            Data[v].append(np.nan)
#        print(v)
#        try :
        #    Data[v].append(float(List[Names.index(v)]))
        #    except IndexError:
        #        print(Names)
        #        print(List)
         #       print(Names.index(v))
         #       print(Data.keys())  
         #       input()
    return Data
    
def parseline(line,delimiter=" "):
    line=line.strip()
    List=line.split(delimiter)   
#    print(List)
    return List

def write_dictionnary(file,keys,List_info,delimiter='\t'):
    """
    Write data into a txt files. Use list file to conserve the ordering of the data 
    INPUT
    -----
    file : str
        output filename
    keys : list
        list of string labeling each set of data (in order)
    List_info : list of np array
        Data to be stored. Each element of the list must have the same length
    OUTPUT
    -----
    None
    
    """
    Dir = os.path.dirname(file)      
    if not os.path.isdir(Dir):
        os.makedirs(Dir)
        
    #dictList is a dictionnary, where each key contain a list of elements of various type
    f = open(file,'w')

    write_header(f,keys)
    
    n=len(List_info[0])
    for i in range(n):
      write_lineDict(f,List_info,i)
    f.close()

def write_a_dict(file,dict_input,delimiter='\t'):
    List_info = []
    keys = dict_input.keys()
    for key in keys:
        List_info.append(dict_input[key])
    write_dictionnary(file,keys,List_info,delimiter=delimiter)
    
       
def write_matrix(file,keys,List_info,dim=2,delimiter='\t'):
    """
    similar to write_dictionnary, but can support multidimensionnal array by flattened it into a 1d array.
    the first dim elements of List_info correspond to the axis labels, the other fields regroups the data. They mus share the same x,y, etc. axis    
    
    """
    #dim =2        
    List = []
    
    #for dim = 2 only for now ! (if greater dimension, might be better to save data in different txt files anyway)
    i = 0
    List.append(np.ndarray.flatten(np.transpose(np.asarray([List_info[i] for j in List_info[i+1]]))))
    List.append(np.ndarray.flatten(np.asarray([List_info[i+1] for j in List_info[i]])))
    for l in List_info[dim:]:
        List.append(np.ndarray.flatten(l))
 
    for l in List:
        print(np.shape(l))
    
    write_dictionnary(file,keys,List,delimiter=delimiter)
 
def write_header(f,keys):
    line=''
    for key in keys:
        line=line+str(key)+'\t'
    line=line+'\n'
    f.write(line)
    
def write_lineDict(f,List,i,delimiter='\t'):
    #List is a List of list that correspond to the data to write
 #   print(List[0])
    line=''
    for L in List:
        line=line+str(L[i])+'\t'
    line=line+'\n'
    f.write(line)

def write_line(f,List,delimiter='\t'):
    #List is a List of list that correspond to the data to write
 #   print(List[0])
    line=''
    for L in List:
        line=line+str(L)+'\t'
    line=line+'\n'
    f.write(line)
#file='/Volumes/labshared3/Stephane/Experiments/Accelerated_grid/2015_08_03/M_2015_08_03/Corr_functions/Corr_spatial_2015_08_03_1_0.txt'
#Header,Data = read_dataFile(file,Hdelimiter='\t',Ddelimiter='\t')
    
    
    
    