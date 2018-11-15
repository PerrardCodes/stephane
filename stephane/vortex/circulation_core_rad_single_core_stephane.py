#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      Joseph
#
# Created:     14/07/2015
# Copyright:   (c) Joseph 2015
# Licence:     <your licence>
#-------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import glob
import os
from matplotlib.colors import Normalize
from scipy import ndimage, interpolate, signal
from scipy import optimize as opt
import cPickle
import random
import scipy.io as sio

class MidpointNormalize(Normalize): # for making sure white is really 0 in the color map
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))

def rm_big(u,v,b=2):
    """ removes spurious data from 2d arrays of x and y velocities by
        converting to length and phaze of the vectors at that point and
        using a threshold filter on these

        INPUTS
        ------
        u : 2D numpy array
        array of x velocities

        v : 2D numpy array
        array of y velocities

        b : int (default 2)
        how many neigbors to consider on either side of the value
        being checked. Note that this many pixels in will NOT be filtered.

        OUTPUT
        ------
        returns two 2D numpy arrays of filtered data, one of x velocities and the
        other of y velocities.
        """
    mag = np.sqrt(u**2+v**2) # going to look at them in terms of magnitude & phaze
    phaze = np.arctan2(v,u)
    ny,nx = u.shape
    for i in xrange(b,ny-b):
        for j in xrange(b,nx-b):
            mag_nbrs = mag[i-b:i+b+1,j-b:j+b+1]
            mag_med = np.median(mag_nbrs)
            x_phaze_nbrs = phaze[i,j-b:j+b+1]
            med_x_phaze = np.median(x_phaze_nbrs)
            y_phaze_nbrs = phaze[i-b:i+b+1,j]
            med_y_phaze = np.median(y_phaze_nbrs)
            if mag[i,j] > mag_med*1.5: # if too long
                mag[i,j] = mag_med
            if np.abs(phaze[i,j]-med_x_phaze) > np.pi/3 or np.abs(phaze[i,j]-med_y_phaze) > np.pi/3: # if too rotated ( more than ~60 deg )
                phaze[i,j] = (med_x_phaze + med_y_phaze)/2.
    ux = mag * np.cos(phaze)
    uy = mag * np.sin(phaze)
    return ux, uy

# next three functions for removing nans from the data, from Stephane

def rm_nans(u):
    """ same as above but only one 2D array """
    i,j = np.where(np.isnan(u))
    for t in zip(i,j):
        u=replace_nan(u,t)
    return u

def replace_nan(u,t,mode='mean'):
    #replace a NaN value by an average of non nan adjacent values (weighted ??)
    # if on a corner , mirror the other side ??
    # u corresponds to the matrix given
    # t corresponds to the indices of nan value to remove in u
    neigh,ind=neighboors(u,t,b=2)
    if mode=='mean':
        u[t]=np.nanmean(neigh)
    if mode=='median':
        u[t]=np.nanmedian(neigh)

    return u
    #ind = where(logical_not(np.isnan(y)))[0]
#y1 = interp(range(len(y)), ind, y[ind])
#y1 = y1[ind[0]:ind[-1]]

def neighboors(u,t,b=1,return_d=False):
    #case b=1. for large b, corner and edge might be increased
    #corner = (t in [(0,0),(ny-1,0),(0,nx-1),(ny-1,nx-1)])
    #edge = ((i in [0,nx-1]) or (j in [0,ny-1])) and (not corner)
    #bulk = not (corner or edge)
    #only d=2 case. might be extended to higher dimension (up to 4 !)
    d=2
    ny,nx=u.shape
    j=t[0]
    i=t[1]
#    neigh=u[j-b:j+b+1,i-b:i+b+1]
        #number of neighbours per line
#    n=2*b+1
#    neigh_vec=np.reshape(neigh,(n**2,)) # 8 neighbours + the central point
#    neigh_vec0=np.concatenate((neigh_vec[:n*b+b],neigh_vec[n*b+b+1:])) # list of the 8 neighbouhrs elements
    boolean=[[(j-p)**2+(i-q)**2<=b*d and (j-p)**2+(i-q)**2>0  for p in range(nx)] for q in range(ny)] # JM: swapped nx and ny from original, outer one is actually y component
    ind = np.where(boolean)
    nan_neigh=np.where(np.isnan(u[ind])) # JM: Changed this to reflect what is probably intended function
 #   print('Number of neighbors : '+str(len(ind[0])))
  #  print('Number of NaN neighbors : '+str(len(nan_neigh[0])))
    #d can be used later to weight the average by closest point
    if return_d:
        d=[(j-tup[0])**2+(i-tup[1])**2 for tup in ind]
        return u[ind],ind,d
    else:
        return u[ind],ind
# above three functions for removing nans from the data, from Stephane

def partial_deriv(plane,direction,dx):
    """
    Finds the partial derivative in a numpy arrray
    INPUT
    -----
    plane : 2D numpy array
    a 2D numpy array giving the values of the function at each point
    direction : string, 'y' or 'x'
    take y or x partial derivative. Assumes positive is up and right
    dx : float
    the spacing between adjascent values in the independent coordinate
    OUTPUT
    ------
    This returns a 2D numpy array of the same dimension as the input array, with
    the value of the desired partial derivative at each point.
    """
    result = np.zeros(plane.shape)
    if direction == 'y':
        for i in xrange(plane.shape[1]): # do each column
            result[0,i] = (plane[0,i] - plane[1,i])/ dx # top edge case. Note indexing is 'upside-down'
            for j in xrange(1,plane.shape[0]-1): # middle cases
                result[j,i] = (plane[j-1,i] - plane[j+1,i])/(2*dx)
            result[plane.shape[0]-1,i] = (plane[plane.shape[0]-2,i] - plane[plane.shape[0]-1,i])/dx # bottom edge case
        return -1*result
    elif direction == 'x':
        for i in xrange(plane.shape[0]): # do each row
            result[i,0] = (plane[i,0] - plane[i,1]) / dx # left edge case
            for j in xrange(1,plane.shape[1]-1): # middle cases
                result[i,j] = (plane[i,j-1] - plane[i,j+1])/(2*dx)
            result[i,plane.shape[1]-1] = (plane[i,plane.shape[1]-2] - plane[i,plane.shape[1]-1]) / dx # right edge case
        return -1*result # -1 because I messed up...
    else:
        print 'Direction not specified or was invalid! please try again using x or y'

def interp_flow_component(x,y,u,piv_roi):
    """
    returns a function where you give x and y coordinates and it gives you interpolated
    values for some component of the flow.

    INPUT
    -----
    x : 1D numpy array
    array of x values

    y : 1D numpy array
    array of y values

    u : 1D numpy array
    array of flow component values

    piv_roi : list
    list describing endpoints of rectangle of the image we really care about.
    Should be in format [x_min,x_max,y_min,y_max]

    OUTPUT
    ------
    Returns the grid of data points, meaning a numpy array of shape (y.len,x.len)
    """
    uindicies = zip(y,x)
    grid_ind_y, grid_ind_x = np.mgrid[piv_roi[2]:piv_roi[3], piv_roi[0]:piv_roi[1]] # gives a grid between specified points (of spacing 1 instead of PIV boxes)
    return interpolate.griddata(uindicies,u,(grid_ind_y, grid_ind_x),method='linear', fill_value=0.) # this interpolates linearly flows onto the grids made above

def circPts(r,phis,cx,cy):
    """takes a radius and phi coordinates and makes cartesian points on a circle

    INPUT
    -----
    r : float
        The radius of your circle

    phis : 1D numpy array
        list of phi values for your circle

    cx : float
        x coordinate of center of the circle

    cy : float
        y coordinate of center of the circle

    OUTPUT
    ------

    returns two 1D numpy arrays, one for x coordinates of circle points and one
    for y coordinates of circle points
    """
    return r*np.cos(phis) + cx, r*np.sin(phis) + cy

def circ_by_rad(r,g,a):
    """ Function to curve fit circulatoin as a function of radial distance from
    a vortex core with a gaussian vorticity profile """
    return g * (1-np.exp(-r**2/a**2))

def ellip_gauss_func(x,y,cx,cy,sx=1.,sy=1.):
    """ an elliptical gaussian bump for convolution, to find the center of the vortex core """
    z = (x-cx)**2/sx**2 + (y-cy)**2/sy**2
    return np.exp(-z)

def gauss_func_cfit((x,y),cx,cy,sx=1.,sy=1.,a=1.):
    """ read somewhere that suggested that we need the first thing as a tuple for curve fit """
    z = (x-cx)**2/sx**2 + (y-cy)**2/sy**2
    return np.exp(-z)
    
def main():
    # needs input
    direc ="/Volumes/labshared/Stephane_lab1/Vortices/2016_01_07/PIV_data/PIVlab_ratio2_W32pix_PIV_vortex_Rjoseph_fps4000_f100mm_tube_d12mm_p30mum_v50/" # NOTE: make sure this is consistent. Some inconsistancy in the file structure/name
    name = direc.split('/')[7] # NOTE: make sure this is consistent. Some inconsistancy in the file structure/name
    frame_begin = 115 # enter the frame range you would like to do the processing on
    frame_end = 150
    frame_by = 1 # How many frames to skip, of form [frame_begin:frame_end:frame_by]
    ppmm = 33.5 # pixels per mm
    radii = np.arange(10,250,15) # radii to calculate the circulation at (in pixels)
    num_pts = 800 # number of points to have on the circles
    position = 'in' # which side of the core is being looked at. in = negative vorticity, out = positive

    # all automatic
    file_dirs = sorted(glob.glob(os.path.join(direc,'*.txt'))) # finding the names of the files data is kept in
    fs = float(name.split('_')[6][3:])
    result_dir = './'+name+'_circ_to_core_DPIV/' # name of the result directory.
    if not os.path.isdir(result_dir):
        os.mkdir(result_dir)
    core_radii = {} # Dictionary which will equate a frame number to the core radius at that frame
    rainbow_cm = cm.gist_rainbow
    rad_list = [] # list of radii found
    if not os.path.isdir(os.path.join(result_dir,'%s_figs/'%position)): # Directory to save the figures in. Is subdir of the result directory
        os.mkdir(os.path.join(result_dir,'%s_figs/'%position))

    # Get a first estimate of the core radius to make a gaussian bump to convolve with the other images
    # Picks out three random frames which are used to get an estimate of the core size.
    rad_est = []
    options = [random.randrange(frame_begin,frame_end+1),random.randrange(frame_begin,frame_end+1),random.randrange(frame_begin,frame_end+1)] # 3 random frames
    for frame_no in options: # need to slice these so we only take frames where it is actually in the frame

        # unload all of the information
        x_array, y_array, ux_array, uy_array = np.loadtxt(file_dirs[frame_no],skiprows=3,usecols=(0,1,2,3),unpack=True,delimiter=',')
        frame = frame_no
        x_array = x_array.astype(int)
        y_array = y_array.astype(int)

        # finding information about which points were taken in PIV
        x_beg = x_array[0]
        x_end = x_array[-1]
        x_by = 16 # NOTE: This is particular. I'm not sure of an easy way to find this when all of them are just in a row. Would probably have to do length of differneces to find a number then multiply
        y_beg = y_array[0]
        y_end = y_array[-1]
        y_by = y_array[1] - y_array[0]
        num_x =(x_end - x_beg)/x_by + 1
        num_y = (y_end - y_beg)/y_by + 1
        piv_roi = [0,x_end,0,y_end]

        # Make the arrays 2D
        plane_y, plane_x = np.mgrid[y_beg:y_end+1:y_by,x_beg:x_end+1:x_by]
        plane_ux = np.transpose(np.reshape(ux_array,(num_x,num_y)))
        plane_uy = np.transpose(np.reshape(uy_array,(num_x,num_y))) # NOTE: x goes from right to left!!!


        # removing NaN and spurious points from the data
        good_plane_ux = rm_nans(plane_ux)
        good_plane_uy = rm_nans(plane_uy)
        best_plane_ux, best_plane_uy = rm_big(good_plane_ux,good_plane_uy)
        best_ux_array = np.transpose(best_plane_ux).flatten()
        best_uy_array = np.transpose(best_plane_uy).flatten()

        # compute curl
        duy_dx = partial_deriv(best_plane_uy,'x',x_by) # NOTE Changed this so no -1 in the dx direction
        dux_dy = partial_deriv(best_plane_ux,'y',y_by)
        plane_curl = duy_dx - dux_dy

        if position == 'out':
            # take maximum of vorticity as center of vortex core
            ext_curl = np.where(plane_curl == plane_curl.max())
            ext_cx = ext_curl[1]
            ext_cy = ext_curl[0]
        else:
            # take minimum of vorticity as center of vortex core
            ext_curl = np.where(plane_curl == plane_curl.min())
            ext_cx = ext_curl[1] # NOTE these are wrt the PIV index NOT the real image pixels
            ext_cy = ext_curl[0]

        # interpolate
        i_ux = interp_flow_component(x_array,y_array,best_ux_array,piv_roi) # returning a grid interpolated to every 1 (real) pixel instead of every 16 which is what we get from PIV
        i_uy = interp_flow_component(x_array,y_array,best_uy_array,piv_roi)
        phis = np.arange(0,2*np.pi,2*np.pi/float(num_pts)) # for plotting circles
        gammas = []
        for rad in radii:
            circx, circy = circPts(rad,phis,ext_cx*x_by,ext_cy*y_by) # circular path to do integral around
            ds = rad*2*np.pi/float(num_pts) # element of arc length
            i_ux_path = ndimage.map_coordinates(i_ux,[circy,circx],order=1) # points of ux interpolated on the path
            i_uy_path = ndimage.map_coordinates(i_uy,[circy,circx],order=1)
            gamma = sum(i_uy_path * np.cos(phis) - i_ux_path * np.sin(phis))*ds # the integral
            gammas.append(gamma)
        gammas = np.asarray(gammas) # found circulations
        real_gammas = gammas * fs * (1/ppmm)**2 # converting to real units, mm^2/s
        real_radii = radii / ppmm # converting to real units, mm
        popt, pcov = opt.curve_fit(circ_by_rad,real_radii,real_gammas,p0=(real_gammas.max(),4.)) # fitting to the correct distribution, give rough first est of core size of 4 mm
        sigma = np.abs(popt[1]) # should be core size in mm
        rad_est.append(sigma) # just need an estimate
        print 'done with estimate %d'%frame_no
    med_rad_est = np.median(rad_est)*ppmm # estimate for mean radius from the 3 random frames
    print 'estimated radius: %.2f mm'%(med_rad_est/ppmm)

    # need to make the thing to convolve
    # making the thing to convolve
    y,x = np.mgrid[0:num_y,0:num_x]
    gauss_bump = ellip_gauss_func(x,y,num_x/2.,num_y/2.,med_rad_est/x_by,med_rad_est/y_by)

    # this is the real processing now
    k = frame_begin
    for file_dir in file_dirs[frame_begin:frame_end+1:frame_by]: # need to slice these so we only take frames where it is actually in the frame

        # unload all of the information
        x_array, y_array, ux_array, uy_array = np.loadtxt(file_dir,skiprows=3,usecols=(0,1,2,3),unpack=True,delimiter=',')
        frame = k
        x_array = x_array.astype(int)
        y_array = y_array.astype(int)

        # finding information about which points were taken in PIV
        x_beg = x_array[0]
        x_end = x_array[-1]
        x_by = 16 # NOTE: This is particular. I'm not sure of an easy way to find this when all of them are just in a row. Would probably have to do length of differneces to find a number then multiply
        y_beg = y_array[0]
        y_end = y_array[-1]
        y_by = y_array[1] - y_array[0]
        num_x =(x_end - x_beg)/x_by + 1
        num_y = (y_end - y_beg)/y_by + 1
        piv_roi = [0,x_end,0,y_end]

        # Make the arrays 2D
        plane_y, plane_x = np.mgrid[y_beg:y_end+1:y_by,x_beg:x_end+1:x_by]
        plane_ux = np.transpose(np.reshape(ux_array,(num_x,num_y)))
        plane_uy = np.transpose(np.reshape(uy_array,(num_x,num_y))) # NOTE: x goes from right to left!!!

        # removing NaN and spurious points from the data
        good_plane_ux = rm_nans(plane_ux)
        good_plane_uy = rm_nans(plane_uy)
        best_plane_ux, best_plane_uy = rm_big(good_plane_ux,good_plane_uy)
        best_ux_array = np.transpose(best_plane_ux).flatten()
        best_uy_array = np.transpose(best_plane_uy).flatten()

        # compute curl
        duy_dx = partial_deriv(best_plane_uy,'x',x_by) # NOTE Changed this so no -1 in the dx direction
        dux_dy = partial_deriv(best_plane_ux,'y',y_by)
        plane_curl = duy_dx - dux_dy

        # convolve with a gaussian bump to find vortex core
        convolution = signal.convolve2d(gauss_bump,plane_curl)

        # find and rescale the convolution center
        if position == 'out':
            conv_ext_indx = np.where(convolution == convolution.max()) # IDK Why but seemed like we had to go to -1 instead of +1 here
            conv_ext_x = (conv_ext_indx[1][0] - num_x-1)*x_by + x_beg
            conv_ext_y = (conv_ext_indx[0][0] - num_y-1)*y_by + y_beg
            # convolution only gives us the positions relative to the center of the image so translate
            ext_cx = conv_ext_x + x_end/2.
            ext_cy = conv_ext_y + y_end/2.
        else:
            conv_ext_indx = np.where(convolution == convolution.min())
            conv_ext_x = (conv_ext_indx[1][0] - num_x-1)*x_by + x_beg
            conv_ext_y = (conv_ext_indx[0][0] - num_y-1)*y_by + y_beg
            # convolution only gives us the positions relative to the center of the image so translate
            ext_cx = conv_ext_x + x_end/2.
            ext_cy = conv_ext_y + y_end/2.

        # sub-pixel interpolation for convolution
        ygrid, xgrid = np.mgrid[0:num_y,0:num_x]
        xgrid_f = xgrid.ravel()*x_by # should be in real pixel location now
        ygrid_f = ygrid.ravel()*y_by
        if position == 'out':
            init_guess = (ext_cx,ext_cy, med_rad_est,med_rad_est,1)
        else:
            init_guess = (ext_cx,ext_cy, med_rad_est,med_rad_est,-1)
        popt, pcov = opt.curve_fit(gauss_func_cfit, (xgrid_f,ygrid_f), plane_curl.ravel(),p0 = init_guess)
        ext_cx_n, ext_cy_n = popt[:2] # should be an even better approximation

        # interpolate
        i_ux = interp_flow_component(x_array,y_array,best_ux_array,piv_roi) # returning a grid interpolated to every 1 (real) pixel instead of every ?? which is what we get from PIV
        i_uy = interp_flow_component(x_array,y_array,best_uy_array,piv_roi)
        phis = np.arange(0,2*np.pi,2*np.pi/float(num_pts)) # for circles
        print 'frame %d'%frame

        pos_dir = '%s_figs/'%position
        gammas = []
        for rad in radii:
            circx, circy = circPts(rad,phis,ext_cx_n,ext_cy_n) # circular path to do integral around
            ds = rad*2*np.pi/float(num_pts) # element of arc length
            i_ux_path = ndimage.map_coordinates(i_ux,[circy,circx],order=1) # points of ux interpolated on the circular path
            i_uy_path = ndimage.map_coordinates(i_uy,[circy,circx],order=1)
            gamma = sum(i_uy_path * np.cos(phis) - i_ux_path * np.sin(phis))*ds # the integral
            gammas.append(gamma)
        gammas = np.asarray(gammas) # found circulations
        real_gammas = gammas * fs * (1/ppmm)**2 # converting to real units, mm^2/s
        real_radii = radii / ppmm # converting to real units, mm
        popt, pcov = opt.curve_fit(circ_by_rad,real_radii,real_gammas,p0=(real_gammas.max(),med_rad_est)) # fitting to the correct distribution
        fig = plt.figure(figsize=(18.0,8.0)) # makes the figure 18x8 inches~
        plt.subplot(221)
        plt.title('ux for frame %d'%frame)
        norm = MidpointNormalize(midpoint=0) # for making sure white is really 0 in the color map
        plt.imshow(best_plane_ux,norm=norm,cmap=cm.bwr,extent=[0,x_end,0,y_end]) # plot of ux

        # plotting circles we took the integrals over
        for rad in radii:
            circx, circy = circPts(rad,phis,ext_cx_n,ext_cy_n)
            circy = y_end - circy # because the origin of the image is at the top but origin for plot is at the bottom
            plt.plot(circx,circy,color=rainbow_cm(float(rad)/radii.max()),zorder=10)
        plt.subplot(222)
        plt.title('uy for frame %d'%frame)
        norm = MidpointNormalize(midpoint=0) # for making sure white is really 0 in the color map
        plt.imshow(best_plane_uy,norm=norm,cmap=cm.bwr,extent=[0,x_end,0,y_end]) # plot of uy

        # plotting circles we took the integrals over
        for rad in radii:
            circx, circy = circPts(rad,phis,ext_cx_n,ext_cy_n)
            circy = y_end - circy # because the origin of the image is at the top but origin for plot is at the bottom
            plt.plot(circx,circy,color=rainbow_cm(float(rad)/radii.max()))

        # plotting calculated circulation as a function of radius of path
        plt.subplot(212)
        fit_dom = np.linspace(0,real_radii.max(),100)
        fit_range = circ_by_rad(fit_dom,popt[0],popt[1])
        plt.plot(fit_dom,fit_range,color='k',zorder=1)
        for i, rad in enumerate(real_radii):
            plt.scatter(rad,real_gammas[i],marker='o',color=rainbow_cm(float(rad)/real_radii.max()),zorder=2)
        plt.xlabel('radius (mm)')
        plt.ylabel('circulation (mm^2/s)')
        plt.xlim(0,real_radii.max()*1.01)
        plt.ylim(real_gammas.min()-np.abs(real_gammas.min())*0.01,real_gammas.max()+np.abs(real_gammas.max())*0.01)
        plt.savefig(os.path.join(result_dir,pos_dir,'circulation_v_rad_frame_%d.png'%frame))
        plt.close()
        sigma = np.abs(popt[1]) # should be core size in mm
        print 'measured core radius (%s): %.2f'%(position,sigma)

        # plotting vorticity with an overlay of the core
        norm = MidpointNormalize(midpoint=0) # for making sure white is really 0 in the color map
        plt.imshow(plane_curl,norm=norm,cmap=cm.bwr,extent=[0,x_end,0,y_end])
        circx, circy = circPts(sigma*ppmm,phis,ext_cx_n,ext_cy_n)
        circy = y_end - circy # because the origin of the image is at the top but origin for plot is at the bottom
        plt.plot(circx,circy,'g')
        plt.scatter([ext_cx_n],[y_end - ext_cy_n],marker='x',color='g')
        plt.title('frame %d'%frame)
        plt.savefig(os.path.join(result_dir,pos_dir,'found_core_overlay_frame_%d.png'%frame))
        plt.close()
        core_radii[frame] = sigma # saving
        rad_list.append(sigma)
        k += frame_by

    # plotting core radius over time and the mean
    t = np.arange(0,len(rad_list))*frame_by
    mean_rad = sum(rad_list)/len(rad_list)
    median_rad = np.median(rad_list)
    top_y = max(rad_list)
    fig = plt.figure(figsize=(18.0,8.0)) # makes the figure 18x8 inches~ (basically able to save it maximized now!)
    plt.plot(t,rad_list,'b.-')
    plt.plot([t[0],t[-1]],[median_rad,median_rad],'r')
    plt.title('%s measured core radius over time'%position)
    plt.xlabel('frame number from beginning of analysis')
    plt.ylabel('core radius (mm)')
    plt.ylim(0,max(rad_list)*1.01)
    plt.xlim(t[0]*0.99,t[-1]*1.01)
    plt.savefig(os.path.join(result_dir,'core_rad_time.png'))
    plt.close()

    # saving the info we collected
    core_rad_file = open(os.path.join(result_dir,'core_radii.p'),'wb')
    cPickle.dump(core_radii,core_rad_file)
    core_rad_file.close()



if __name__ == '__main__':
    main()




