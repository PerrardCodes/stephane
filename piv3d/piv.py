# -*- coding: utf-8 -*-
"""
Created on Fri Aug 25 12:18:54 2017

@author: danjr
"""

import pims
import numpy as np
import openpiv_SP.process as process
import openpiv_SP.validation as validation
import openpiv.filters
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import pickle
import pandas as pd

import os
#import comps
from skimage import exposure
import fluids2d.backlight as backlight

import stephane.display.graphes as graphes


from mpl_toolkits.axes_grid1 import make_axes_locatable


class PIVDataProcessing:
    '''
    Class for PIV computations for a cine.
    '''

    def __init__(self,parent_folder,cine_name,window_size=32,overlap=16,search_area_size=64,frame_diff=1,name_for_save=None,maskers=None,crop_lims=None,dx=1,dt_orig=1):

        # Metadata
        self.cine_filepath = parent_folder+cine_name+'.cine'
        self.parent_folder=parent_folder
        self.cine_name=cine_name
        #self.num_frames_orig = len(pims.open(self.cine_filepath))
        self.maskers=maskers
        self.crop_lims = crop_lims
        self.cine_frame_shape=None

        # PIV Parameters
        self.window_size=window_size
        self.overlap=overlap
        self.search_area_size=search_area_size
        self.frame_diff = frame_diff # pairs of frames separated by how many frames

        # Scaling, for reference -- results are still stored in [pixels / frame] !
        self.dx=dx
        self.dt_orig=dt_orig # "original" dt - between frames in the cine
        self.origin_pos=None

        self.dt_ab = self.dt_orig*self.frame_diff # dt between frames a and b

        # For saving the results
        self.flow_field_res_filepath = None # will store a path to the numpy array with the flow field results
        if name_for_save is None:
            self.name_for_save = self.cine_name
        else:
            self.name_for_save = name_for_save

    def save(self):
        self.ff = None # So the actual flowfield isn't stored in this pickled object
        pickle.dump(self,open(self.parent_folder+self.name_for_save+'.pkl','wb'))

    def process_frame(self,frame_a,frame_b,s2n_thresh=1.3):
        frame_a = frame_a.astype(np.int32)
        frame_b = frame_b.astype(np.int32)

        u,v,sig2noise = process.extended_search_area_piv( frame_a, frame_b, window_size=self.window_size, overlap=self.overlap, dt=1, search_area_size=self.search_area_size,sig2noise_method='peak2peak' )
        u, v, mask = validation.sig2noise_val( u, v, sig2noise, threshold = s2n_thresh )

        return u,v

    def run_analysis(self,a_frames=None,save=True,s2n_thresh=1.3,bg_n_frames=None):

        c=pims.open(self.cine_filepath)
        if a_frames is None:
            a_frames = np.arange(0,len(c)-self.frame_diff)

        self.cine_frame_shape = c.frame_shape

        im0 = c.get_frame(0)
        (Nx,Ny)=im0.shape
        Nt = (a_frames[-1]-a_frames[0])//self.frame_diff

        Vol = np.zeros((Nt,Nx,self.frame_diff))
        N = Nt*self.frame_diff
        j0=600
        for i in range(N):
            Vol[i//self.frame_diff,:,i%self.frame_diff]=c.get_frame(i)[:,j0]
        cnew = [Vol[i,:,:int(self.frame_diff/2)] for i in range(Nt)]

        self.cine_frame_shape = cnew[0].shape
        print('shape image (x,z) plane : '+str(cnew[0].shape))
        a_frames = np.arange(0,Nt-1)

        if bg_n_frames is not None:
            bg = backlight.construct_bg_image(c,n_frames=bg_n_frames,usemax=False)
            plt.figure()
            plt.imshow(bg)


        '''
        Define function to crop the image, based on the limits in crop_lims
        '''
        crop_lims = self.crop_lims
        if crop_lims is not None:
#            def crop(im):
#                im[0:crop_lims[2],:] = 0
#                im[crop_lims[3]:-1,:] = 0
#                im[:,0:crop_lims[0]] = 0
#                im[:,crop_lims[1]:-1] = 0
#
#                plt.figure()
#                plt.imshow(im)
#                plt.show()
#
#                return im
            def crop(im): return im[crop_lims[0]:crop_lims[1],crop_lims[2]:crop_lims[3]]
        else:
            def crop(im): return im

        #flow_field = init_flowfield(len(a_frames),np.shape(crop(c[0])),self.window_size,self.overlap)
        flow_field = init_flowfield(len(a_frames),np.shape(crop(cnew[0])),self.window_size,self.overlap)

        '''
        Store some processing parameters
        '''
        #self.window_coordinates = process.get_coordinates(np.shape(crop(c[0])),self.window_size,self.overlap)
        self.window_coordinates = process.get_coordinates(np.shape(crop(cnew[0])),self.window_size,self.overlap)
        self.s2n_thresh=s2n_thresh # threshold used in PIV analysis
        self.a_frames=a_frames # which original frames are the "a" frames
        self.dt_ab = self.dt_orig*self.frame_diff # dt between frames a and b
        self.dt_frames = self.dt_orig * np.median(np.diff(self.a_frames)) # ASSUMES THAT ALL A FRAMES ARE EVENLY SPACED!
        try:
            self.a_frame_times = [c.frame_time_stamps[i] for i in a_frames] # time of the a frames
        except:
            print('frame time stamps not available')
            self.a_frame_times = a_frames*self.dt_frames

        for aii,ai in enumerate(a_frames):
            print('file '+str(self.name_for_save)+', frame a: '+str(ai))

            # get the two frames
            #frame_a = crop(c[ai].astype(np.int32))
            #frame_b = crop(c[ai+self.frame_diff].astype(np.int32))

            frame_a = cnew[ai]#c[ai]
            frame_b = cnew[ai+1]#c[ai+self.frame_diff]

#            frame_a = exposure.adjust_gamma(frame_a.astype(float),1.3).astype(np.int32)
#            frame_b = exposure.adjust_gamma(frame_b.astype(float),1.3).astype(np.int32)

            if bg_n_frames is not None:
                frame_a = frame_a-bg
                frame_b = frame_b-bg

            frame_a = frame_a.astype(np.int32)
            frame_b = frame_b.astype(np.int32)

            if self.maskers is not None:
                for masker in self.maskers:
                    frame_a = masker.mask_frame(ai,frame_a).astype(np.int32)
                    frame_b = masker.mask_frame(ai,frame_b).astype(np.int32)

            frame_a = crop(frame_a)
            frame_b = crop(frame_b)

            # run the PIV analysis
            u,v,sig2noise = process.extended_search_area_piv( frame_a, frame_b, window_size=self.window_size, overlap=self.overlap, dt=1, search_area_size=self.search_area_size,sig2noise_method='peak2peak' )
            u, v, mask = validation.sig2noise_val( u, v, sig2noise, threshold = s2n_thresh )

            # store the velocity fields
            flow_field[aii,:,:,0] = u
            flow_field[aii,:,:,1] = v

            if (save==True)&(ai%10==0): # save every some frames, just in case
                print(flow_field.shape)
                np.save(self.parent_folder+self.name_for_save+'_flowfield.npy',flow_field)
                self.flow_field_res_filepath = self.parent_folder+self.name_for_save+'_flowfield.npy'
                self.save()

        if save==True:
            np.save(self.parent_folder+self.name_for_save+'_flowfield.npy',flow_field)
            self.flow_field_res_filepath = self.parent_folder+self.name_for_save+'_flowfield.npy'
            self.save()

        return flow_field

    def load_flowfield(self):
        '''
        Return an instance of the PIVDataAnalysis class containing the SCALED
        dataset.
        '''
        flow_field = PIVDataAnalysis(load_scaledflowfield(self.parent_folder+self.name_for_save+'_flowfield.npy',self.dx,self.dt_ab))
        print('dx : ')
        print(self.dx)
        print('dt_ab: ')
        print(self.dt_ab)
        return flow_field

    def associate_flowfield(self):
        '''
        So the flowfield can be accessed through the .data attribute.

        If this object is going to be saved, the .data attribute will first be
        cleared so as to not save the actual flowfield with the metadata
        contained in this object.
        '''
        self.data = self.load_flowfield()

#    def indx2time(self,idx):
#        c=pims.open(self.cine_filepath)
#        time = [c]


class PIVDataAnalysis:
    '''
    Class to work with processed PIV data. Data is ALREADY SCALED when
    initialized!

    The main attribute is "ff" (flowfield), a 4D numpy array:
        local_inst_velocity_component = ff[frame,row,column,direction]
    '''
    def __init__(self,ff):
        '''
        ff here is ALREADY SCALED given dx and dt!
        '''
        self.ff = ff

    def get_frame(self,i):
        '''
        Return one frame from the flowfield array (4D -> 3D)
        '''
        return self.ff[i,:,:,:]

    def show_frame(self,i,bg='speed',ax=None):
        frame = self.get_frame(i)
        ax=show_frame(frame[:,:,0],frame[:,:,1],bg=bg,ax=ax)
        return ax

    def average_frames(self,frames):
        return np.nanmean(self.ff[frames,:,:,:],axis=0)

class TurbulenceParams:
    '''
    Class for mapping forcing conditions to turbulence parameters
    '''

    def __init__(self,metadata,data):
        '''
        metadata is a dataframe indexed with the dict data. it has the forcing
        conditions associated with each entry in data
        '''

###############################################################################

'''
FUNCTIONS FOR DATA I/O
'''

def read_control_csv(csv_filepath):
    '''
    Read the .csv file that controls which cases are to be analyzed.
    '''

    all_cases = pd.read_csv(csv_filepath)
    all_cases = all_cases.loc[all_cases['use_now']==1]

    amplitudes = all_cases['A'].unique()
    freqs = all_cases['freq'].unique()

def combine_linear_data(data,geometries,dy=0.001):
    '''
    List of numpy arrays and list of corresponding geometries. The arrays will
    be stacked vertically with some interpolation.
    '''

    from scipy.interpolate import interp1d

    '''
    First find the vertical extent of all the geometries.
    '''
    min_y = np.nan
    max_y = np.nan
    for g in geometries:
        min_y = np.nanmin([min_y,np.min(g.y)])
        max_y = np.nanmax([max_y,np.max(g.y)])

    print(min_y)
    print(max_y)
    new_y = np.arange(min_y,max_y,dy)
    new_f = np.nan * new_y

    for di in np.arange(len(data)):
        d = data[di]
        g = geometries[di]

        print(np.shape(d))
        print(np.shape(g.y))

        min_y = np.min(g.y)
        max_y = np.max(g.y)

        interpolated_y = new_y[(new_y>=min_y)&(new_y<=max_y)]
        interpolator = interp1d(g.y,d)
        interpolated_f = interpolator(interpolated_y)

        #new_f[(new_y>=min_y)&(new_y<=max_y)] = np.nanmean([new_f[(new_y>=min_y)&(new_y<=max_y)],interpolated_f])
        new_f[(new_y>=min_y)&(new_y<=max_y)] = interpolated_f

    return new_y,new_f

def rotate_data_90(ff,n=1):
    '''
    Rotate the images in a 4d matrix 90 deg ccw.

    Can't use the 'axes' parameter in numpy.rot90 since that version of numpy
    is not compatible with trackpy.
    '''

    for _ in range(n):

        # first move the time dimension out of the way
        ff = np.swapaxes(ff,0,1)
        ff = np.swapaxes(ff,1,2)

        # rotate about the plane of the first two axes
        ff = np.rot90(ff)

        # put back the time axis at the front
        ff = np.swapaxes(ff,2,1)
        ff = np.swapaxes(ff,1,0)

        # flip the velocity components
        ff = ff[:,:,:,[1,0]]

        # correct the direction of the u velocity
        ff[:,:,:,0] = -1*ff[:,:,:,0]

    return ff

def init_flowfield(num_frames,frame_shape,window_size,overlap):
    '''
    Initialize the 4d numpy array to store the flow field
    '''
    field_shape = process.get_field_shape(frame_shape,window_size,overlap)
    flow_field = np.zeros([num_frames,field_shape[0],field_shape[1],2]) # [frame, row, column, velocity component]
    return flow_field

def scale_flowfield(ff,dx,dt_ab):
    '''
    Convert [pixels/frame] to [m/s]
    '''
    return ff * dx / dt_ab

def load_scaledflowfield(filepath,dx,dt):
    '''
    Load the flowfield stored at filepath, and scale it given dx and dt.

    Returns a numpy array.
    '''
    ff=np.load(filepath)
    return scale_flowfield(ff,dx,dt)

def load_processed_PIV_data(parent_folder,name):
    '''
    Load the pickled metadata object and call the .associate_flowfield method
    so the stored flowfield is in the .ff attribute.
    '''
    p = pickle.load(open(parent_folder+name+'.pkl','rb'))
    p.parent_folder = parent_folder
    p.associate_flowfield() # the numpy data can now be accessed with p.data.ff[...]
    return p

def clip_flowfield(ff,thresh):
    ff2 = ff.copy()
    #ff2[ff>thresh] = thresh
    #ff2[ff<-1*thresh] = -1*thresh
    ff2[ff>thresh] = np.nan
    ff2[ff<-1*thresh] = np.nan
    return ff2

def fill_nans(x):
    '''
    Use pandas to forward fill then backfill nan values in a series
    '''
    return pd.Series(x).fillna(method='ffill').fillna(method='bfill').values

def fill_nans_3d(x):
    '''
    Forward/backward fill nans in a 3d array along the first axis
    '''
    num_rows = np.shape(x)[1]
    num_cols = np.shape(x)[2]

    for yi in range(num_rows):
        for xi in range(num_cols):
            x[:,yi,xi] = fill_nans(x[:,yi,xi])

    return x

def fill_nans_nd(x):
    s = np.shape(x)

    to_it = s[1:]
    idxs = [range(i) for i in to_it]

    for i in zip(idxs):
        print(i)
        print(np.shape(x[:,i]))
        x[:,i] = fill_nans(x[:,i])

def check_concatenation(p):

    first_job_name = p.cine_name+'_job0'
    p_0 = pickle.load(open(p.parent_folder+first_job_name+'.pkl'))
    p_0.parent_folder = p.parent_folder
    p_0.name_for_save = first_job_name
    p_0.associate_flowfield()

    is_good = p.data.ff[0,0,0,0]==p_0.data.ff[0,0,0,0]

    return is_good



###############################################################################

'''
FUNCTIONS FOR FLOWFIELD ANALYSIS
'''

def compute_frame_shear(u,v):
    uxx = np.gradient(np.gradient(u,axis=1),axis=1)
    vyy = np.gradient(np.gradient(v,axis=0),axis=0)
    shear = 0.5 * ( uxx**2 + vyy**2 )
    return shear

def compute_gradients(ff,dx=1):
    '''
    Given a (4D) 3D flowfield, return a (5D) 4D gradient field indexed by
    [(time,)row,column,velocity_dir,gradient_dir]

    for example, du/dx is found in ff_gradients[...,0,1]
    '''

    ff_gradients = np.zeros(np.shape(ff)+(2,))

    for i in [0,1]:
        for j in [0,1]:
            ff_gradients[...,i,j] = np.gradient(ff[...,i],axis=j+1) / dx

    return ff_gradients

def compute_turbulent_fluctuations(ff):
    mean_flow = np.nanmean(ff,axis=0)
    return ff - mean_flow

def compute_3d_outofplane_velocity(all_data,g,dz):
    '''
    Use the continuity equation for incompressible flow to compute the w
    component of velocity for a scanning PIV dataset. all_data contains the 3d,
    2-component velocity, indexed as
        vel = all_data[time,y_pos,x_pos,piv_plane,vel_component]
    where y_pos and x_pos are the positions in the original 2d images,
    piv_plane gives the location of the plane in the perpendicular direction,
    and vel_component in (0,1) corresponds to the u or v velocity in-plane with
    the original images.

    g should be the GeometryScaler class for a given 2d PIV plane (it's taken
    to be identical for all the planes), and dz is the spacing between planes
    (also taken to be constant).

    It returns an array with_new_data, which is identical to all_data except
    there are 3 components of velocity (ie vel_component can be in (0,1,2)).

    Since the continuity equation only gives dw/dz, after integration through
    the planes from 0, each "line" of w values perpendicular to the original
    PIV planes is adjusted by some value such that the average value is 0; this
    is obviously not suitable for some cases.
    '''

    # compute the known gradients of velocity
    dudx = np.gradient(all_data[...,0],axis=2) / g.dx
    dvdy = np.gradient(all_data[...,1],axis=1) / g.dx

    # get the third one from incompressibility
    dwdz = 1 - dudx - dvdy

    # integrate dw/dz to get w, shifted by some arbitrary amount
    w = np.cumsum(dwdz,axis=3) * dz

    # normalize so that over each line through the planes, the mean out-of-plane
    # velocity is 0
    w = w - np.nanmean(w,axis=(0,1,2))

    # make a new array that contains this w component
    with_new_data = np.append(all_data,np.moveaxis(np.array([w]),0,-1),axis=4)

    return with_new_data


def basic_field_calcs(ff):
    '''
    Return a dict of common fields (mean flow, fluctuations, etc)
    '''

    results = {}

    results['mean_flow'] =np.nanmean(ff,axis=0)
    results['fluc'] = ff-mean_flow
    results['u_rms'] = np.sqrt( np.nanmean( (results['fluc'][:,:,:,0])**2,axis=0) + np.nanmean( (results['fluc'][:,:,:,1])**2,axis=0) )
    results['inst_speed'] = np.linalg.norm(ff,ord=2,axis=3)
    results['meanflow_speed'] = np.sqrt((results['mean_flow'][:,:,0])**2 + (results['mean_flow'][:,:,1])**2)
    results['dudx'] = np.gradient(results['fluc'][:,:,:,0],axis=1)
    results['dvdy'] = np.gradient(results['fluc'][:,:,:,1],axis=2)

    return results

def dissipation_direct(ff,g,nu=1e-6):

    dudx = np.gradient(ff[...,0],axis=2)/g.dx
    dudy = np.gradient(ff[...,0],axis=1)/g.dx
    dvdx = np.gradient(ff[...,1],axis=2)/g.dx
    dvdy = np.gradient(ff[...,1],axis=2)/g.dx

    epsilon_direct = 4*nu* ( np.nanmean(dudx**2,axis=0) + np.nanmean(dvdy**2,axis=0) + np.nanmean(dudx*dvdy,axis=0) + 0.75*np.nanmean((dvdx+dudy)**2,axis=0) )

    return epsilon_direct


def compute_radial_and_angular_velocity(vel,g):
    '''
    Convert to polar coordinates and find the radial and angular velocities,
    about (0,0) as define in g.
    '''
    R = np.sqrt(g.X**2+g.Y**2)
    Theta = np.arctan2(g.Y,g.X)
    U_r = vel[...,0]*np.cos(Theta) + vel[...,1]*np.sin(Theta)
    U_theta = -1*vel[...,0]*np.sin(Theta) + vel[...,1]*np.cos(Theta)
    return U_r,U_theta,R,Theta

def compute_DLL_and_DTT(U_r,U_theta,R,Theta,g,n_theta=200,n_r=51,r_max=0.075):
    '''
    Given u_r and u_theta, compute the quantities used to find the longitudinal
    and transverse structure functions.

    DLL and DTT should be averaged over axes 0 and 2 to get just a function of
    r.
    '''

    # how many points in time
    n_t = np.shape(U_r)[0]

    theta_vec = np.linspace(0,2*np.pi,n_theta)
    r_half_vec = np.linspace(0,r_max,n_r)
    r_vec = r_half_vec*2

    DLL = np.zeros([n_t,n_r,n_theta])
    DTT=np.zeros_like(DLL)
    #fig = plt.figure()
    #ax = fig.add_subplot(111)
    for thi,theta in enumerate(theta_vec):
        '''
        Go from theta=0 to theta=2pi
        '''
        for ri,r in enumerate(r_half_vec):

            xi = np.argmin(np.abs(g.x-r*np.cos(theta)))
            yi = np.argmin(np.abs(g.y-r*np.sin(theta)))

            xi_neg = np.argmin(np.abs(g.x+r*np.cos(theta)))
            yi_neg = np.argmin(np.abs(g.y+r*np.sin(theta)))

            #ax.plot([xi,xi_neg],[yi,yi_neg])

            DLL[:,ri,thi] = (U_r[:,yi,xi] + U_r[:,yi_neg,xi_neg])**2 # (minus a negative)
            DTT[:,ri,thi] = (U_theta[:,yi,xi] + U_theta[:,yi_neg,xi_neg])**2

    #plt.show()

    return DLL,DTT,r_vec,theta_vec

def compute_DLL_and_DTT_fromCenter(U_r,U_theta,R,Theta,g,n_theta=200,n_r=51,r_max=0.05,sep_center=0.2):
    '''
    Given u_r and u_theta, compute the quantities used to find the longitudinal
    and transverse structure functions.

    DLL and DTT should be averaged over axes 0 and 2 to get just a function of
    r.
    '''

    # how many points in time
    n_t = np.shape(U_r)[0]

    theta_vec = np.linspace(0,2*np.pi,n_theta)
    r_vec = np.linspace(0,r_max,n_r)

    DLL = np.zeros([n_t,n_r,n_theta])
    DTT=np.zeros_like(DLL)
    #fig = plt.figure()
    #ax = fig.add_subplot(111)
    for thi,theta in enumerate(theta_vec):
        '''
        Go from theta=0 to theta=2pi
        '''
        for ri,r in enumerate(r_vec):

            xi = np.argmin(np.abs(g.x-r*np.cos(theta)))
            yi = np.argmin(np.abs(g.y-r*np.sin(theta)))

            xi_neg = np.argmin(np.abs(g.x+r*np.cos(theta)*sep_center))
            yi_neg = np.argmin(np.abs(g.y+r*np.sin(theta)*sep_center))

            #ax.plot([xi,xi_neg],[yi,yi_neg])

            DLL[:,ri,thi] = (U_r[:,yi,xi] + U_r[:,yi_neg,xi_neg])**2 # (minus a negative)
            DTT[:,ri,thi] = (U_theta[:,yi,xi] + U_theta[:,yi_neg,xi_neg])**2

    #plt.show()

    return DLL,DTT,r_vec*(1.+sep_center),theta_vec

def compute_autocorr(vel,R,Theta,g,n_theta=200,n_r=51,r_max=0.075,sep_center=0.5):
    '''
    Given u_r and u_theta, compute the quantities used to find the longitudinal
    and transverse structure functions.

    DLL and DTT should be averaged over axes 0 and 2 to get just a function of
    r.
    '''

    # how many points in time
    n_t = np.shape(vel)[0]

    theta_vec = np.linspace(0,2*np.pi,n_theta)
    r_vec = np.linspace(0,r_max,n_r)

    B_LL = np.zeros([n_t,n_r,n_theta])
    #fig = plt.figure()
    #ax = fig.add_subplot(111)
    for thi,theta in enumerate(theta_vec):
        '''
        Go from theta=0 to theta=2pi
        '''
        for ri,r in enumerate(r_vec):

            xi = np.argmin(np.abs(g.x-r*np.cos(theta)))
            yi = np.argmin(np.abs(g.y-r*np.sin(theta)))

            xi_neg = np.argmin(np.abs(g.x+r*np.cos(theta)*sep_center))
            yi_neg = np.argmin(np.abs(g.y+r*np.sin(theta)*sep_center))

            #ax.plot([xi,xi_neg],[yi,yi_neg])

            disp_norm = [np.cos(theta),np.sin(theta)]
            if r==0:
                disp_norm=np.array([1,1])/np.sqrt(2)

            vel_center_comp = np.sum(np.tile(disp_norm,(n_t,1))*vel[:,yi_neg,xi_neg,:],axis=1)
            vel_comp = np.sum(np.tile(disp_norm,(n_t,1))*vel[:,yi,xi,:],axis=1)

            #prod =

            B_LL[:,ri,thi] = vel_center_comp * vel_comp

    #plt.show()

    return B_LL,r_vec*(1+sep_center),theta_vec


def compute_autocorr_fromCenter(vel,R,Theta,g,n_theta=200,n_r=51,r_max=0.05,):
    '''
    Given u_r and u_theta, compute the quantities used to find the longitudinal
    and transverse structure functions.

    DLL and DTT should be averaged over axes 0 and 2 to get just a function of
    r.
    '''

    # how many points in time
    n_t = np.shape(vel)[0]

    theta_vec = np.linspace(0,2*np.pi,n_theta)
    r_vec = np.linspace(0,r_max,n_r)

    B_LL = np.zeros([n_t,n_r,n_theta])
    #fig = plt.figure()
    #ax = fig.add_subplot(111)
    for thi,theta in enumerate(theta_vec):
        '''
        Go from theta=0 to theta=2pi
        '''
        for ri,r in enumerate(r_vec):

            xi = np.argmin(np.abs(g.x-r*np.cos(theta)))
            yi = np.argmin(np.abs(g.y-r*np.sin(theta)))

            xi_mid = np.argmin(np.abs(g.x))
            yi_mid = np.argmin(np.abs(g.y))

            #ax.plot([xi,xi_neg],[yi,yi_neg])

            # x before y here to match the velocity which is [u,v]
#            disp = np.array([g.x[xi]-g.x[xi_mid],g.y[yi]-g.y[yi_mid]])
#            norm = np.linalg.norm(disp)
#            if norm==0:
#                disp_norm = np.array([1.,1.])
#            else:
#                disp_norm = disp/norm
            disp_norm = [np.cos(theta),np.sin(theta)]
            if r==0:
                disp_norm=np.array([1,1])/np.sqrt(2)

            vel_center_comp = np.sum(np.tile(disp_norm,(n_t,1))*vel[:,yi_mid,xi_mid,:],axis=1)
            vel_comp = np.sum(np.tile(disp_norm,(n_t,1))*vel[:,yi,xi,:],axis=1)

            #prod =

            B_LL[:,ri,thi] = vel_center_comp * vel_comp
            #B_TT[:,ri,thi] = (U_theta[:,yi,xi]*U_theta[:,yi_neg,xi_neg]*-1)

    #plt.show()

    return B_LL,r_vec,theta_vec

###############################################################################

'''
FUNCTIONS FOR VISUALIZATION
'''

def add_fieldimg_to_ax(ff_reduced,g,ax,time=None,slice_dir=None,vel_dir=None,cmap_other=None,vmin=-0.2,vmax=0.2):
    '''
    Add a 2d image to an axes with the appropriate scaling.

    Use for either a snapshot or time average, or spatially-reduced composite
    image with appropriate slice_dir defined
    '''

    if vel_dir=='horizontal':
        cmap = 'PuOr'
    elif vel_dir=='vertical':
        cmap = 'seismic'
    elif vel_dir==None:
        cmap = 'viridis'

    if cmap_other is not None:
        cmap = cmap_other

    if slice_dir=='vertical':
        cb=ax.imshow(ff_reduced.transpose(),origin='upper',vmin=vmin,vmax=vmax,aspect='auto',cmap=cmap,extent=[0,max(time),g.im_extent[2],g.im_extent[3]])
    elif slice_dir=='horizontal':
        cb=ax.imshow(ff_reduced,vmin=vmin,vmax=vmax,aspect='auto',cmap=cmap,extent=[g.im_extent[0],g.im_extent[1],0,max(time)],origin='lower')
    elif slice_dir==None:
        cb=ax.imshow(ff_reduced,vmin=vmin,vmax=vmax,cmap=cmap,extent=[g.im_extent[0],g.im_extent[1],g.im_extent[2],g.im_extent[3]],origin='upper')

    return cb

def composite_image_plots(ff,g,time,row_lims,col_lims,vmin=-0.2,vmax=0.2):

    '''
    Velocity along center column
    '''

    fig_column = plt.figure(figsize=(16,9))
    ax1=fig_column.add_subplot(1,2,1)
    ax2=fig_column.add_subplot(1,2,2)

    #ff = ff - np.nanmean(ff,axis=0)

    add_fieldimg_to_ax(np.nanmean(ff[:,:,col_lims[0]:col_lims[1]+1,0],axis=2),g,ax1,time=time,slice_dir='vertical',vel_dir='horizontal',vmin=vmin,vmax=vmax)
    add_fieldimg_to_ax(np.nanmean(ff[:,:,col_lims[0]:col_lims[1]+1,1],axis=2),g,ax2,time=time,slice_dir='vertical',vel_dir='vertical',vmin=vmin,vmax=vmax)

    [a.set_xlabel('Time [s]') for a in [ax1,ax2]]
    ax1.set_ylabel('Position along vertical column [m]')
    ax1.set_title('Horizontal velocity')
    ax2.set_title('Vertical velocity')

    '''
    Velocity along center span
    '''

    fig_row = plt.figure(figsize=(16,9))
    ax1=fig_row.add_subplot(1,2,1)
    ax2=fig_row.add_subplot(1,2,2)

    add_fieldimg_to_ax(np.nanmean(ff[:,row_lims[0]:row_lims[1]+1,:,0],axis=1),g,ax1,time=time,slice_dir='horizontal',vel_dir='horizontal',vmin=vmin,vmax=vmax)
    add_fieldimg_to_ax(np.nanmean(ff[:,row_lims[0]:row_lims[1]+1,:,1],axis=1),g,ax2,time=time,slice_dir='horizontal',vel_dir='vertical',vmin=vmin,vmax=vmax)
    [a.set_xlabel('Position along horizontal span [m]') for a in [ax1,ax2]]
    ax1.set_ylabel('Time [s]')
    ax1.set_title('Horizontal velocity')
    ax2.set_title('Vertical velocity')

    return fig_column,fig_row

def vertical_percentile_distributions(ff,g,time,row_lims,col_lims,other_percentiles=[10,25,75,90]):

    fig=plt.figure(figsize=(9,4))
    ax1=fig.add_subplot(131)
    ax2=fig.add_subplot(132,sharey=ax1,sharex=ax1)
    ax3=fig.add_subplot(133,sharey=ax1)

    [a.axvline(0,color='gray',alpha=0.3) for a in [ax1,ax2,ax3]]


    ax1.plot(np.nanmean(np.nanpercentile(ff[:,:,col_lims[0]:col_lims[1],0],50,axis=0),axis=1),g.y,color='k',alpha=1)
    for other_percentile in other_percentiles:
        ax1.plot(np.nanmean(np.nanpercentile(ff[:,:,col_lims[0]:col_lims[1],0],other_percentile,axis=0),axis=1),g.y,color='k',alpha=0.5)

    ax2.plot(np.nanmean(np.nanpercentile(ff[:,:,col_lims[0]:col_lims[1],1],50,axis=0),axis=1),g.y,color='k',alpha=1)
    for other_percentile in other_percentiles:
        ax2.plot(np.nanmean(np.nanpercentile(ff[:,:,col_lims[0]:col_lims[1],1],other_percentile,axis=0),axis=1),g.y,color='k',alpha=0.5)

    mean_flow=np.nanmean(ff,axis=0)
    fluc = ff-mean_flow
    u_rms = np.sqrt( np.nanmean( (fluc[:,:,:,0])**2,axis=0) + np.nanmean( (fluc[:,:,:,1])**2,axis=0) )

    ax3.plot(np.nanmean(np.linalg.norm(mean_flow[:,col_lims[0]:col_lims[1],:],ord=2,axis=2),axis=1),g.y,color='k',alpha=0.5,label='''$\sqrt{\overline{U}^2 + \overline{V}^2}$''')
    ax3.plot(np.nanmean(u_rms[:,col_lims[0]:col_lims[1]],axis=1),g.y,color='k',alpha=1,label='''$\sqrt{\overline{(u-\overline{U})^2}+\overline{(v-\overline{V})^2}}$''')

    ax1.set_title('u distribution')
    ax2.set_title('v distribution, mean is'+str(np.nanmean(np.nanmean(np.nanpercentile(ff[:,:,col_lims[0]:col_lims[1],1],50,axis=0),axis=1))))
    ax3.set_title('fluctuations')
    ax3.legend()
    ax1.set_ylabel('vertical position [m]')
    [a.set_xlabel('velocity [m/s]') for a in [ax1,ax2,ax3]]

    return fig



def plot_both_components(vel,time=None):
    fig=plt.figure(figsize=(9,5))
    ax1=fig.add_subplot(121)
    ax2=fig.add_subplot(122)

    mean_flow = np.nanmean(vel,axis=0)
    print(np.shape(mean_flow))

    if time is None:
        time = np.arange(0,np.shape(vel)[0])

    ax1.axhline(0,color='k',alpha=0.2)
    ax1.plot(time,vel[:,0],'--',label='u',color='gray',lw=0.8)
    ax1.plot(time,vel[:,1],'-.',label='v',color='gray',lw=0.8)
    ax1.plot(time,np.sqrt((vel[:,0]-mean_flow[0])**2+(vel[:,1]-mean_flow[1])**2),color='r',lw=1,label='''$\sqrt{ u'^2 + v'^2}$''')
    ax1.plot(time,np.sqrt((vel[:,0])**2+(vel[:,1])**2),color='k',lw=2,label='''$\sqrt{ u^2 + v^2}$''')
    ax1.legend()

    ax2.axvline(0,color='k',alpha=0.2)
    ax2.axhline(0,color='k',alpha=0.2)
    ax2.scatter(vel[:,0],vel[:,1],c=time,alpha=0.3)
    ax2.set_ylabel('v')
    ax2.set_ylabel('u')

def overlay_vectors_frame(p,frame):

    fig = plt.figure()
    ax=fig.add_subplot(111)

    indx_in_cine = p.a_frames[frame]

    cine=pims.open(p.cine_filepath)
    ax.imshow(cine[indx_in_cine],cmap='gray')

    '''
    Show the mask
    '''
    mask = p.maskers[0].create_mask(indx_in_cine)
    ax.imshow(mask,alpha=0.5)


def show_mean_flows(ff):
    fig=plt.figure(figsize=(8,8))
    ax1=fig.add_subplot(2,2,1)
    ax2=fig.add_subplot(2,2,2)
    ax3=fig.add_subplot(2,2,3)
    ax4=fig.add_subplot(2,2,4)

    ax1.imshow(np.nanmean(ff[:,:,:,0],axis=0))
    ax2.imshow(np.nanmean(ff[:,:,:,1],axis=0))
    ax3.imshow(np.nanstd(ff[:,:,:,0],axis=0))
    ax4.imshow(np.nanstd(ff[:,:,:,1],axis=0))

def show_frame(u,v,bg='speed',ax=None):
    if ax is None:
        fig=plt.figure(); ax=fig.add_subplot(111)

    if bg=='speed':
        bg = np.sqrt(u**2+v**2)
    elif bg=='shear':
        bg = compute_frame_shear(u,v)

    ax.matshow(bg,vmin=0)
    ax.quiver(u,v,color='white',alpha=0.5)
    plt.show()
    plt.pause(0.1)
    return ax

def convert_ax_to_cbar(ax,cmap,vmin,vmax):

    import matplotlib as mpl
    pos = ax.get_position().bounds
    ax.set_position([pos[0],pos[1]+pos[3]/2,pos[2],(pos[3])/5])

    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    cb1 = mpl.colorbar.ColorbarBase(ax, cmap=cmap,
                                    norm=norm,
                                    orientation='horizontal')
    return cb1,ax



class PIVComparisonsFigures:
    '''
    Class to create the figures needed for each comparison of PIV cases
    '''

    def __init__(self,nrows,ncols,figsize=(11,8),max_speed=1.0,vmin=-.5,vmax=0.5,epsilon_max=1000,legend_axi=None):

        self.nrows = nrows
        self.ncols = ncols
        self.figsize = figsize
        self.legend_axi = legend_axi
        self.left_axs = np.arange(0,nrows*ncols,ncols)
        self.bottom_axs = np.arange(ncols*(nrows-1),ncols*nrows,1)

        self.max_speed = max_speed
        self.vmin = vmin
        self.vmax = vmax
        self.epsilon_max = epsilon_max

        self.fig_names = ['mean','u_rms','intns','U','V','isotropy','epsilon']

        self.figs = {}
        self.axs = {}
        for name in self.fig_names:
            self.figs[name],self.axs[name] = plt.subplots(self.nrows,self.ncols,figsize=self.figsize)
            self.axs[name] = self.axs[name].flatten()

    def add_case(self,ai,ff,g,time):

        mean_flow=np.nanmean(ff,axis=0)
        fluc = ff-mean_flow
        u_rms = np.sqrt( np.nanmean( (fluc[:,:,:,0])**2,axis=0) + np.nanmean( (fluc[:,:,:,1])**2,axis=0) )
        meanflow_speed = np.sqrt((mean_flow[:,:,0])**2 + (mean_flow[:,:,1])**2)
        isotropy_rms = np.sqrt( np.nanmean( (fluc[:,:,:,0])**2,axis=0)) / np.sqrt( np.nanmean( (fluc[:,:,:,1])**2,axis=0))
        epsilon = dissipation_direct(fluc,g)

        # Add the images to the figures
        add_fieldimg_to_ax(np.log10(u_rms / meanflow_speed),g,self.axs['intns'][ai],time=time,slice_dir=None,vel_dir=None,vmin=-1,vmax=1,cmap_other='coolwarm')
        add_fieldimg_to_ax(meanflow_speed,g,self.axs['mean'][ai],time=time,slice_dir=None,vel_dir=None,vmin=0,vmax=self.max_speed)
        add_fieldimg_to_ax(u_rms,g,self.axs['u_rms'][ai],time=time,slice_dir=None,vel_dir=None,vmin=0,vmax=self.max_speed)
        add_fieldimg_to_ax(mean_flow[:,:,0],g,self.axs['U'][ai],vel_dir='horizontal',vmin=self.vmin,vmax=self.vmax)
        add_fieldimg_to_ax(mean_flow[:,:,1],g,self.axs['V'][ai],vel_dir='vertical',vmin=self.vmin,vmax=self.vmax)
        add_fieldimg_to_ax(np.log2(isotropy_rms),g,self.axs['isotropy'][ai],slice_dir=None,vel_dir=None,vmin=-2,vmax=2,cmap_other='PRGn')
        add_fieldimg_to_ax(epsilon*100**2,g,self.axs['epsilon'][ai],slice_dir=None,vel_dir=None,vmin=0,vmax=self.epsilon_max,cmap_other='viridis')

        # update the limits
        [g.set_axes_limits(a[ai]) for a in list(self.axs.values())]

        # clean up the axis labels etc
        if ai not in self.left_axs:
            pass
           # [a[ai].yaxis.set_ticklabels([]) for a in list(self.axs.values())]
        else:
            [a[ai].set_ylabel('y [m]') for a in list(self.axs.values())]

        if ai not in self.bottom_axs:
            pass
            #[a[ai].xaxis.set_ticklabels([]) for a in list(self.axs.values())]
        else:
            [a[ai].set_xlabel('x [m]') for a in list(self.axs.values())]

    def update_limits(self,ai,xlims,ylims):
        [a[ai].set_xlim(xlims) for a in list(self.axs.values())]
        [a[ai].set_ylim(ylims) for a in list(self.axs.values())]

    def add_rect(self,ai,x,y,width,height,color,ls,lw=2):
        for a in list(self.axs.values()):
            r = matplotlib.patches.Rectangle([x,y],width,height,edgecolor=color,ls=ls,lw=lw,fill=False)
            a[ai].add_patch(r)

    def add_text(self,ai,x,y,text):
        [a[ai].text(x,y,text,ha='right',va='top',color='k',fontsize=8) for a in list(self.axs.values())]

    def tight_layout(self):
        [f.tight_layout() for f in list(self.figs.values())]

    def add_legends(self):

        # intensity
        cbar,ax_cbar = convert_ax_to_cbar(self.axs['intns'][self.legend_axi],'coolwarm',-1,1)
        cbar.set_ticks([-1,0,1])
        cbar.set_ticklabels([0.1,1,10])
        ax_cbar.set_title('turbulence intensity')

        # mean flow speed
        cbar,ax_cbar = convert_ax_to_cbar(self.axs['mean'][self.legend_axi],'viridis',0,self.max_speed)
        ax_cbar.set_title('mean flow speed [m/s]')

        # u_rms
        cbar,ax_cbar = convert_ax_to_cbar(self.axs['u_rms'][self.legend_axi],'viridis',0,self.max_speed)
        ax_cbar.set_title('$u_\mathrm{rms}$ [m/s]')

        # mean U
        cbar,ax_cbar = convert_ax_to_cbar(self.axs['U'][self.legend_axi],'PuOr',self.vmin,self.vmax)
        ax_cbar.set_title('$\overline{U}$ [m/s]')

        # mean V
        cbar,ax_cbar = convert_ax_to_cbar(self.axs['V'][self.legend_axi],'seismic',self.vmin,self.vmax)
        ax_cbar.set_title('$\overline{V}$ [m/s]')

        # isotropy
        cbar,ax_cbar = convert_ax_to_cbar(self.axs['isotropy'][self.legend_axi],'PRGn',-2,2)
        cbar.set_ticks([-2,0,2])
        cbar.set_ticklabels([.25,1,4])
        ax_cbar.set_title('''$\sqrt{ \overline{u'^2} } / \sqrt{ \overline{w'^2} }$''')

        # dissipation
        cbar,ax_cbar = convert_ax_to_cbar(self.axs['epsilon'][self.legend_axi],'viridis',0,self.epsilon_max)
        ax_cbar.set_title('$\epsilon$ [cm^2/s^2]')

    def remove_axes(self,ai):
        [a[ai].set_axis_off() for a in list(self.axs.values())]

    def save_figs(self,figfolder,prefix):
        for key in list(self.figs.keys()):
            fig = self.figs[key]
            fig.savefig(figfolder+prefix+'_'+key+'.pdf')

def comparison_figures_from_df(meta,nrows,ncols,figsize=(11,8),max_speed=1.0,vmin=-.5,vmax=0.5,legend_axi=None):


    #Figs = piv.PIVComparisonsFigures(3,10,figsize=(17,10),max_speed=0.2,vmin=-.15,vmax=0.15,legend_axi=0)
    Figs = PIVComparisonsFigures(nrows,ncols,figsize=(11,8),max_speed=max_speed,vmin=vmin,vmax=vmax,legend_axi=legend_axi)
    skip_ax=[]


    fig_fft = plt.figure()
    ax_fft = fig_fft.add_subplot(111)

    '''
    Loop through each case and add the plots to the figures
    '''

    vmin = -0.2
    vmax = 0.2

    rect_size = 0.005

    colors=['r','g','b','cyan',]

    C_dict = {}
    lags_dict = {}
    for i in meta.index:

        parent_folder = meta.loc[i,'parent_folder']
        case_name = meta.loc[i,'case_name']
        offset = meta.loc[i,'offset']
        color = colors[i]

        ai = i
        for si in [Figs.legend_axi] + skip_ax:
            if ai>=si:
                ai = ai+1

        #label = 'T = '+str(meta.loc[i,'period'])+' ms\n$\phi$ = '+str(int(100*meta.loc[i,'on_portion']))+'/100'
        #if meta.loc[i,'description'] is not None:
        #    label = label+'\n'+meta.loc[i,'description']

        p = pickle.load(open(parent_folder+case_name+'.pkl'))
        p.parent_folder = parent_folder
        p.name_for_save = case_name
        p.associate_flowfield()

        if meta.loc[i,'need2rotate']:
            p.data.ff=piv.rotate_data_90(p.data.ff)
            p.data.ff=piv.rotate_data_90(p.data.ff)
            p.data.ff=piv.rotate_data_90(p.data.ff)
        ff = p.data.ff

        g_orig = fluids2d.geometry.GeometryScaler(dx=p.dx,im_shape=(1,1),origin_pos=offset,origin_units='pix')
        g = fluids2d.geometry.create_piv_scaler(p,g_orig)

        print(np.shape(ff))
        print(g.im_shape)
        print(g.im_extent)

        time = np.arange(0,np.shape(ff)[0]) * p.dt_frames

        lim = 0.03
        center_rows,center_cols = g.get_coords(np.array([[rect_size,-rect_size],[-rect_size,rect_size]]))

        '''
        Filter the velocity field
        '''
        ff=piv.clip_flowfield(ff,5)

        '''
        Mean flow and fluctuations
        '''




        '''
        Add to the figures
        '''
        Figs.add_case(ai,ff,g,time)
        #Figs.update_limits(ai,[-0.1,0.05],[-0.12,0.1])
        #Figs.update_limits(ai,[-0.10,0.10],[-0.11,0.1])
        #Figs.add_text(ai,0.,.09,label)
        #Figs.add_rect(ai,-rect_size,-rect_size,rect_size*2,rect_size*2,color=color)


    Figs.tight_layout()
    Figs.add_legends()
    [Figs.remove_axes(a) for a in skip_ax]
    return Figs

def gen_artificial_polar_flow(n_t=10000,n_rows=41,n_cols=51):
    '''
    Artificial polar field to test the code

    Flow is proportional to sqrt(X^2+Y^2) plus some noise
    '''
    nrows = 41
    ncols = 51
    X,Y = np.meshgrid(np.arange(ncols),np.arange(nrows))
    X = X-ncols/2
    Y = Y-nrows/2
    ff = np.sqrt(X**2+Y**2)
    ff = (np.moveaxis(fluc*np.ones((1000,2,nrows,ncols)),1,-1)+1) * np.random.rand(1000,nrows,ncols,2)
    #ff = ff - np.nanmean(fluc,axis=0)

    return ff


if __name__ == '__main__':

    '''
    Run the module directly to do a single PIV case.
    '''

    '''
    Setup the processing
    '''

    # which case to process
    parent_folder = '/media/stephane/'
    date = '20181010'


#    cine_name = 'OS/Documents and Settings/Stephane/Documents/Data_buffer/20181010/PIV3dscan_nikon50mm_f1kHz_A800mV_offsetm2800mV_rotator'

    cine_name = 'OS/Documents and Settings/Stephane/Documents/Data_buffer/20181010/PIV3dscan_nikon50mm_f1kHz_A800mV_offsetm2800mV_4pumpsOn'

    #parent_folder = r'\\Mae-deike-lab3\c\Users\Luc Deike\data_comp3_C\180726\\'
    #parent_folder = r'F:\from_comp3_D\171203\\'
    #parent_folder = comps.cf('comp3d')+r'180803\\'


    # dx and dt for the cine file
    #dx = 3.3240951317E-05 # viewA
    #dx = 8.4191870904E-05 # viewB
    dx =  0.35E-03#7.3469357156E-05
    dt_orig = 1./40000

    crop_lims=None
    pre_constructed_masker = None
    t0 = 0
    a_frames = np.arange(t0,40000,1)
    frame_diff = 40

    window_size = (32,20)
    overlap =  (16,8)
    search_area_size=window_size


    data_folder = 'DATA/Experimental_data/Turbulence3d/'+date+'/test3d3c/'
    save_folder = data_folder+ 'window'+str(window_size[0])+'x'+str(window_size[1])+'/'

    if not os.path.isdir(parent_folder+save_folder):
        os.makedirs(parent_folder+save_folder)

#    os.mkdir(parent_folder+save_folder)

    '''
    Run the PIV processing
    '''
    processing = PIVDataProcessing(parent_folder,cine_name,name_for_save=save_folder,dx=dx,dt_orig=dt_orig,frame_diff=frame_diff,crop_lims=crop_lims,maskers=None,window_size=window_size,overlap=overlap,search_area_size=search_area_size)

#Compute flow field
    processing.run_analysis(a_frames=a_frames,save=True,s2n_thresh=1.2,bg_n_frames=None)

# Load existing flowfield
    processing.associate_flowfield()

    ff = processing.data.ff
#    ff = clip_flowfield(ff,4.)

    #convert 2d to 3d data
    (Nt,Nx,Ny,Nc) = ff.shape
    ff = np.reshape(ff,(int(Nt/frame_diff),frame_diff,Nx,Ny,Nc))
    print(ff.shape)
    print(np.sum(np.isnan(ff.flatten()))/np.prod(ff.shape))

    #compute mean_flow
    mean_flow = np.nanmean(ff,axis=0)
    mean_flow_speed = np.linalg.norm(mean_flow,axis=2)
    mean_speed = np.nanmean( np.sqrt(ff[...,0]**2 + ff[...,1]**2 ), axis=0)
    fluc = ff - mean_flow
    u_rms = np.sqrt(np.nanmean( fluc[...,0]**2+fluc[...,1]**2 ,axis=0) )

    d = len(mean_flow.shape)-1
    '''┴
    Show the mean flow
    '''
    (Nz,Nx,Ny,Nc) = mean_flow.shape

    Nz = Nz-1
    dz = 45/15.5
    z0 = -6*dz
    z = np.arange(-Nz/2,Nz/2)*dz-z0
    x = np.arange((Nx-1)/2,-(Nx-1)/2-1,-1)*processing.dx*1E3*overlap[1]+205
    y = np.arange(-(Ny-1)/2,(Ny-1)/2+1)*processing.dx*1E3*overlap[0]

    [Y,X] = np.meshgrid(y,x)
    #representation in the (x,y) plane

    #definition of mask
    x0 = 175
    y0 = 130
    width = 120
    height = 50

    count=0
    print(X.shape)
    print(ff[0,0,:,:,0].shape)

    fig,axs=plt.subplots(1,2,figsize=(6,5),sharey=True); axs=axs.flatten()
    cax=[]
    for i in range(2):
        divider = make_axes_locatable(axs[i])
        cax.append(divider.append_axes('right', size='15%', pad=0.1))


    for j in range(2,18):
        c=[]
        for (i,data) in zip(range(2),[mean_speed[j,:,:],u_rms[j,:,:]]):
            axs[i].clear()
            cb=axs[i].pcolormesh(X,Y,data,vmin=0,vmax=0.6)
#            cb=axs[i].pcolormesh(X,Y,ff[0,j,:,:,1],vmin=0,vmax=0.3)
            c.append(fig.colorbar(cb, cax=cax[i], orientation='vertical'))
            axs[i].set_xlabel('$x$ (mm)')

            p=patches.Rectangle((x0,y0),width,height,facecolor='w',edgecolor='r')
            axs[i].add_patch(p)

        axs[0].set_title(r'$y$ = '+str(int(z[j]))+' mm')
        #axs[1].set_title(r'$\bar u_{rms}$')
        axs[0].set_ylabel('$z$ (mm)')
        c[0].set_label(r'$\bar u $ (m/s)')
        c[1].set_label(r'$\bar u_{rms}$ (m/s)')

        plt.pause(0.1)

        filename = parent_folder+save_folder+'zdimension/'+'Fluc_'+str(count)
    	#if not os.path.isdir(parent_folder+save_folder):
        #	os.makedirs(parent_folder+save_folder)
        graphes.save_fig(1,filename,frmt='png',dpi=300,overwrite=False)
        plt.pause(0.1)
        count+=1


    #representation in the plane (z,y)
    [X,Z] = np.meshgrid(x,z)
    print(Z.shape)

    if d==2:
    	fig,axs=plt.subplots(1,2,figsize=(12,8)); axs=axs.flatten()
    	cmaps = ['PuOr','seismic']
    	for i in [0,1]:
            cb=axs[i].imshow(mean_flow[:,:,i],vmin=-.7,vmax=0.7,cmap=cmaps[i])
    if d==3:
        (Nz,Nx,Ny,Nc) = mean_flow.shape
        cmaps = ['PuOr','seismic']
        count=0
        fig,axs=plt.subplots(1,2,figsize=(16,6)); axs=axs.flatten()
        for j in range(Ny-5,0,-1):

            c = []
            for (i,data,v) in zip([0,1],[mean_flow[:,:,j,0],mean_speed[...,j]],[-0.4,-0.4]):
                divider = make_axes_locatable(axs[i])
                cax = divider.append_axes('right', size='10%', pad=0.1)

                cb=axs[i].pcolormesh(X,Z,data,vmin=v,vmax=0.4)

#                graphes.legende('x (mm)','y (mm)','Mean flow')

                c.append(fig.colorbar(cb, cax=cax, orientation='vertical'))

            c[0].set_label('<uz> (m/s)')
            c[1].set_label('<u> (m/s)')


            axs[0].set_title('Mean flow, z = '+str(y[j])+' mm')
            axs[1].set_title('Mean speed')

            axs[0].set_ylabel('y (mm)')
            axs[0].set_xlabel('x (mm)')
            axs[1].set_xlabel('x (mm)')

            #plt.show()
            plt.pause(0.1)

            filename = parent_folder+save_folder+'norm_'+str(count)
            graphes.save_fig(1,filename,frmt='png',dpi=300,overwrite=False)
            plt.pause(0.1)
            count+=1


    '''
    Show the turbulence fluctuations
    '''
    fig = plt.figure()
    ax=fig.add_subplot(111)
    cb=ax.imshow(u_rms,vmin=0,vmax=1)
    fig.colorbar(cb)

    '''
    Show the turbulence intensity
    '''
    fig=plt.figure()
    ax=fig.add_subplot(111)
    cb=ax.imshow(np.log10(u_rms/mean_speed),vmin=-1,vmax=1,cmap='seismic')
    fig.colorbar(cb)
#
    '''
    Show the isotropy
    '''
    isot = np.sqrt( np.nanmean( fluc[...,0]**2,axis=0) ) / np.sqrt( np.nanmean( fluc[...,1]**2,axis=0) )
    fig = plt.figure()
    ax = fig.add_subplot(111)
    cb = ax.imshow(np.log10(isot),vmin=-1,vmax=1,cmap='PRGn')
    fig.colorbar(cb)

    '''
    Histogram
    '''
    fig = plt.figure()
    ax = fig.add_subplot(111)
    bins = np.linspace(-1,1,501)
    for d in [0,1]:
        vel = ff[...,d].flatten()
        vel = vel[~np.isnan(vel)]
        hist,_ = np.histogram(vel,bins=bins)
        ax.plot(bins[1:],hist)
    ax.axvline(processing.dx/processing.dt_ab,color='k')
    ax.axvline(-1*processing.dx/processing.dt_ab,color='k')

#
#    '''
#    Make an animation
#    '''
#    fig = plt.figure()
#    ax = fig.add_subplot(111)
#
#    dvdx = np.gradient(ff[...,1],axis=2)
#    dudy = np.gradient(ff[...,0],axis=1)
#    vort = dvdx-dudy
#    for fi in range(len(a_frames)):
#        ax.clear()
#        ax.imshow(vort[fi,:,:])
#        plt.show()
#        plt.pause(0.5)
