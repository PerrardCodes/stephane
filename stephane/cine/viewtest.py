from numpy import *
import v4d
from scipy import ndimage

#Make something 4D that looks cool
x, y, z = ogrid[:128, :128, :128]

dat = []

levels = [-0.8, 0.1, 0.6]
thickness = 0.1
scale = 0.05

for phase in  2 * pi * arange(0, 1, 1./16):
#    s = sin(x*scale - phase) * sin((y)*scale + phase)  * sin(z*scale - phase) + 0.5*sin((x+y+z) * scale)
    s = sin(x*scale) * sin((y)*scale)  * sin(z*scale) + 0.5*sin((x+y+z) * scale + phase)
    X = zeros(s.shape) 
    
    for l in levels: X += 0.05 * (s > l - thickness/2.) * (s < l + thickness/2.)
    dat.append(ndimage.gaussian_filter(X, 0.5))  #Smooth it a bit.


#Show it
v4d.show(array(dat)) #Data clipped to (0, 1) by views4d if it is float/double type
