

import stephane.display.graphes
import stephane.mdata.M_manip as M_manip
import stephane.display.vfield as vfield
date = '2015_12_28'

#load the first data of the day
M = M_manip.load_Mdata_serie(date,1,0)

nx,ny,nt = M.shape()
Dirname = './Velocity_module_'+M.Id.get_id() #not used now
vfield.make_2dmovie(M,0,nt/2,fignum=1)