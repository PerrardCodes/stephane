# -*- coding: utf-8 -*-
"""
Created on Tue May 26 09:57:15 2015
Generate all python structures to store the data associated to one particular date
Sdata : associated to each cine file. Contains the time stamp, spatial scale, origin in space, origin in time
Mdata : for each Dir of measurements, after PIV processing for instance
@author: stephane
"""

import matplotlib
matplotlib.use('Agg')

import stephane.mdata.Sdata_manip as Sdata_manip
import stephane.mdata.M_manip as M_manip
import stephane.tools.browse as browse

import argparse
import os.path

#script_dir = os.path.dirname(os.path.realpath(__file__))
parser = argparse.ArgumentParser(description="Create Sdata and Mdata objects associated to cine files and PIV measurements")
parser.add_argument('-d',dest='date',default='1987_03_19',type=str,help='date to be processed. Python will look for it in all the folders specified in file_architecture.py')
parser.add_argument('-f',dest='folder',default=None,type=str,help='base folder to be processed. Python will look for cine files inside this specified folder')
parser.add_argument('-s',dest='start',default=None,type=int, help='start processing index of the cinefile List')
parser.add_argument('-e',dest='end',default=None,type=int, help='end processing index of the cinefile List')
args = parser.parse_args()


#import argparse
#global M
#script_dir = os.path.dirname(os.path.realpath(__file__))
#parser = argparse.ArgumentParser(description="Generates Sdata and Mdata for the specified date, load them in the calling environnement")
#parser.add_argument(dest='date',type=str, help='Date of data to be processed, string format')
#args = parser.parse_args()

def load(date):
    """
    Import data from a given date
    INPUT
    -----
    date : str
    OUTPUT
    -----
    Mlist : list of Mdata 
        class Mdata in stephane.mdata.Mdata
        inherited classes : Mdata_PIVlab, Mdata_pyPIV
    """
   # dataDir='/Volumes/labshared3/Stephane/Experiments/Accelerated_grid/2015_03_24/PIV_bv_hp_X25mm_Zm30mm_fps20000_H1000mm_zoom_zoom_S100mm/PIVLab_logSample_16pix/'    
   # Sdata_manip.Sdata_gen_day(date)
  #  Slist=Sdata_manip.load_Sdata_day(date)      
    Mlist=M_manip.load_Mdata_serie(date)
  # dataDir='/Volumes/labshared3/Stephane/Experiments/Accelerated_grid/2015_03_21/PIV_bv_hp_Zm30mm_fps10000_H1000mm_zoom_S100mm/PIVlab_logSample/'
   # M=M_manip.M_gen(Slist[17],dataDir,1,'PIVlab')
  #  dataDir='/Volumes/labshared3/Stephane/Experiments/Accelerated_grid/2015_03_21/PIV_bv_hp_Zm30mm_fps1000_H1000mm_zoom_S100mm/PIVlab_16pix/'
  #  M=M_manip.M_gen(Slist[18],dataDir,1,'PIVlab')

    return Mlist

def process(date):
    setattr(args,'date',date)
    main()
    
def process_list(date_list):
    for date in date_list:
        process(date)
    
def main():
    Sdata_manip.Sdata_gen_day(args.date)
  
    M_manip.Measure_gen_day(args.date)    
  
  #  date = args.date
  #  Mlist = M_manip.load_Mdata_serie(date,data=False)
  #  for M in Mlist:
  #      M.write(data=False)
        
#    locals()['Slist_'+date] = Sdata_manip.load_Sdata_day(date)   
#    locals()['Mlist_'+date] = data_load(date)
        
if __name__ == '__main__':
    main()
