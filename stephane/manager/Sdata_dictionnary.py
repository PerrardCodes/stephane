# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 16:55:10 2015

@author: stephane
"""

#set a dictionnary with all the existing experiments
#this dictionnary should contain the following list of parameters :

#        date
#        self.index=index
#        self.mindex=mindex
#        self.type="Accelerated grid"
#        self.who="SPerrard"

#     self.Xplane=browse.get_number(file,"_X","mm")
#        self.Zplane=browse.get_number(file,"_Z","mm")
#        self.fps=browse.get_number(file,"_fps","_")
#        self.H0=browse.get_number(file,"_H","mm")
#        self.stroke=browse.get_number(file,"_S","mm")
        
 #       self.typeview='sv'
  #      self.typeplane='hp'
        
        # bounds ?
        
    #type of data : bubbles or PIV measurement, 2d or 3d ? do it in literal way
    # for each, the set of measurement that has been done on it ?
        
    #comments ??
        
    #location of the cinefile
    #location of the Datafile ?
    #has been processed or not ?
import glob
import os

import stephane.manager.file_architecture as file_architecture
    
import stephane.tools.rw_data as rw_data
import stephane.tools.browse as browse

import stephane.mdata.Sdata_manip as Sdata_manip



key_dict=['date','index','H0','fps','stroke','Xplane','Zplane','typeview','typeplane','dataType','fileCine','processed','type','who']
default={'dataType':'2D_PIVmovie','processed':'No'}
#default={'dataType':'2D_Bubblemovie','processed':'No'}
rootDir='/Volumes/labshared3/Stephane/Experiments/Accelerated_grid'
#rootDir='/Volumes/labshared/Stephane_lab1'
#parcours l'arborescence de dossier pour indexer tous les films de manips existants ??
#Yes !!!

def generate_dictionnary():
    #for all the subdir in rootDir, find all the cine files, generate a dictionnary.
    failure = []
    for rootDir in file_architecture.Dir_data:
        print(rootDir)
        
        dataDirList,n=browse.get_dirList(rootDir[:-1],'20*',True)
         
        print(dataDirList)
        general_dict=[]
        start=True
        for Dir in dataDirList:
            print(Dir)
            date=browse.get_end(Dir,rootDir)
    
            print(date)
            if True:
            #try:
                dirdata = file_architecture.get_dir(date)
                paramList = glob.glob(dirdata+'/*param.txt')
                print(dirdata+'/*param.txt',paramList)               
                for filename in paramList:
                    os.remove(filename)
#                input()
                
                Sdata_manip.Sdata_gen_day(date)
                
                dictList=dictionnary(date)
            
                if not dictList==[]:
                    write(dictList,Dir,date)
    
                    if start:
                        start=False
                        general_dict=dictList
                    else:
                        for i,elem in enumerate(dictList):
                        #concatene the existing general dictionnnary with the current dictionnary
                            general_dict[i]=general_dict[i] + elem
            else:
#            except:
                failure = failure + [date]
                print(date + 'failed')
                
        #generate a root dictionnary from any cine files info (!)
        date = 'to_2016_10_11'
        print('')
        print('')
        print('**********************')
        print(general_dict)
        write(general_dict,rootDir,date)
    
    print('Cinefiles that fail to load :' + str(failure))


def write(dictList,Dir,date):
    
    if len(dictList)>0:
        SaveDir=Dir+'/Dictionnary'
        print(SaveDir)
        if not os.path.isdir(SaveDir):
            os.makedirs(SaveDir)
        filename=SaveDir+'/Experiment_overview_'+date+'.txt'
        print(dictList)
        
        rw_data.write_dictionnary(filename,key_dict,dictList)
    
    

def update():
    pass
    #update the dictionnary

def dictionnary(date):
    #from the indicating date, generate a dictionnary from all the present cine file
    Slist=Sdata_manip.load_Sdata_day(date)
    print("++++++++++++++++++++++")
    print(Slist)
    if Slist==[] or Slist is None:
        return []
    else:
        dictList=Slist_to_dict(Slist)    
        return dictList

def update_dictionnary(date):
    #Add the missing data to an existing dictionnary : thus checked if the 
    return None
    
def Slist_to_dict(Slist):
    List_info=[[] for key in key_dict]
    for S in Slist:
        List_info=Sdata_to_list(S,List_info)
        
    return List_info
        
def Sdata_to_list(S,List_info):
    #Slist is a list of Sdata : should be converted to a dictionnary, with the keys given in ref_dict :
    for key in key_dict:
        value=find_attr(S,key)
        if value=='':
            print(key+'attribute is missing ! White space instead')
        List_info[key_dict.index(key)].append(value)

    return List_info  

def find_attr(S,name):
    if name=='dataType':
        return find_dataType(S)
    if name=='processed':
        return processed(S)
    #Look, by order of preference, in : Sdata, Sdata.id, Sdata.param
    lsearch=[S.param, S.Id, S]
    l=lsearch.pop()
    while (not lsearch==[]) and (not hasattr(l,name)):
        l=lsearch.pop()
    
    if hasattr(l,name):
        return getattr(l,name)
    else:
        #look in the default register
        if name in default.keys():
            return default[name]
        else:
            print("S and is subobject have no attribute "+name+' ,'+S.fileCine)
            print('Arbitrary set to Unknown')
#            val=input('Manual input for "'+name+'"')
            return 'Unknown'
            
#Auxiliary functions due to the lack of information in the Sdata object
def find_dataType(S):
    s=browse.get_string(S.fileCine,'/PIV','.cine')
    if not s=='':
        return '2D_PIVmovie'

    s=browse.get_string(S.fileCine,'Bubbles','.cine')
    if not s=='':
        return '2D_bubbles'
    
    return 'Unknown'
    
def processed(S):
    dataDirList,n=browse.get_dataDir(S.fileCine)   
    if n==0:
        return 'No'
    else:
        return 'Yes'#dataDirList
    
def main():
    generate_dictionnary()
#    date='2015_03_21'
#    dictionnary(date)

if __name__ == '__main__':
    main()