

import stephane.manager.Data_representation as Data_representation
import stephane.display.graphes as graphes
import stephane.mdata.Sdata_manip as Sdata_manip

#read the Sdata, ie only the headers of the experiments
Slist = Sdata_manip.load_all()

#plot them graphycally. One graph correspond to one experimental configuration
#return a dict of the figures, with an associated default filename base on the x and y legends
figs=Data_representation.graphic(Slist)
    
#save the figures in a subfolder
graphes.save_figs(figs,savedir='./Figures/',frmt='png')
