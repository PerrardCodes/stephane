# -*- coding: utf-8 -*-
"""
Created on Tue Sep 29 13:12:55 2015

@author: stephane
"""

# catalog of useful tools to characterize a turbulent flow. Some are not yet implemented. 
#Should be upgraded on the way

import stephane.analysis.decay as decay
import stephane.analysis.Fourier as Fourier
import stephane.analysis.early_time as early_time
import stephane.analysis.corr as corr
import stephane.analysis.taylor_scale as taylor_scale

import pylab as plt
import os.path
import stephane.tools.browse as browse

"""
Main programm to run the different kind of data analysis :
    Call the following modules :
    decay : decay of the total energy in the field
    corr : correlation functions in time
        compute_Ct : autocorrelation function in time
        spatial_correlation : correlation in space, by taking N random pair of indices in each slot of distance d,d+Delta between them
    Fourier : fourier transform in 1d or 2d in space
        display_fft_vs_t : compute the 1d Fourier spectrum over time, and fit by a K-spectrum. Eventually plot the Energy contained in the turbulence cascade as a function of time
    early_time : temporal evolution of the flow (isotropy, homogeneity)
        isotropy : plot the isotropy coefficient between the two axis as a function of time
    taylor_scale : Calculation of the Taylor microscale, ie the typical length scale of spatial variation of the strain tensor 
        taylor_scale : currently based on a spatial averaged of U∂_{x,y}U.
        to be implemented : from the curvature of the spatial correlation function for d -> 0 (!)
        
    Old measurement technics, to be updated :
        Distribution of velocity or its increments for various spatial scales
        Decaying energy for each vertical line
        Vertical transport of turbulent kinetic energy over time
            
    To be added :
        Computation of the mean flow from a set of same experimental parameters runs ("smart" alignement needed to adjust slightly different camera positions)
        Vorticity calculation, distribution of vorticity
"""

def measure(Mlist,start):
    functions=[decay.decay,Fourier.display_fft_vs_t,early_time.isotropy,corr.compute_Ct,corr.spatial_correlation]#]#]#c]#corr.compute_Ct]#]#,,,,taylor_scale.taylor_scale]            
    for i,fun in enumerate(functions):
        print(str(fun.__name__))
        iterate(Mlist,fun,start)
       
#    input()
#    iterate(Mlist,corr.spatial_correlation,2)
    #starter method to launch the principal functions used to characterize a turbulent flow    
    #fun : name of the function to be launched
        #decay : compute the Total energy as a function of time. Use some smoothing in time
           # input()
#        display_fft_vs_t(M,'1d',Dt=50,label=label)

        #corr module : contain two different one, correlation function in time and correlation functions in space
            # Time correlation. Compute time correlation functions, and the associated characteristic correlation time tc over time
            # Spatial correlation. Compute spatial correlation functions and save it in txt file for a later use. Return the characteristic correlation length lc over 

        #Fourier : contain measurements of both 1d and 2d spatial spectrum at any given time. Smoothing option in time
        #it is really relevant for measurements from the bottom view

def iterate(Mlist,fun,i):
    exclude= ['2015_09_22_2_0','2015_09_22_0_2','2015_09_22_5_0','2015_09_22_3_2']

    count=0
    for j,M in enumerate(Mlist):
        if M.id.get_id() not in exclude:
            fignum = chose_figure(M,i)
            label = chose_label(M)
            count+=1
        #chose the label as a function of         
            if fignum>0:
                fun(M,label=label,display=True,fignum=fignum)
                print(j)                
                print(M.id.get_id())                

def chose_figure(M,i=0):
    fx_bound = 0.115
#    fx_max = [0.2,0.115]
#    fx_min = [0.115,0.07]

    if M.param.fx>fx_bound:
        fignum=1+i*2
    else:
        fignum=1+i*2
        
    return fignum
    
def chose_label(M):
    base=os.path.basename(M.dataDir)
  #  i=browse.get_string(base,'olymer')
    print(base)
    c=browse.get_string(base,start='_',end='ppm',from_end=True)
    print(c)
    if len(c)>0:
        print(c)
        if int(c)==100:
            color = 'b'     
            symbol = ''
        if int(c)==200:
            color = 'r'     
            symbol = ''
        label = color + symbol
        
    else:
        label = ''
    return label
        
#plt.close("all")

#measure(Mlist_2015_06_04[1:],1)
            
#measure([Mlist_2015_09_28[0]],2)
#measure([Mlist_2015_09_17[0]]+Mlist_2015_09_17[2:5],1)

#measure(Mlist_2015_08_14,1)

#measure(Mlist_2015_09_22,1)

#measure(Mlist_2015_09_25,3)
#
#measure(Mlist_2015_09_28,1)
def main():
    measure(Mlist_2015_11_06,1)

if __name__ == '__main__':
    main()