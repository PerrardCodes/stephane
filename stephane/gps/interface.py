import numpy as np
import pylab as plt
import glob
import os

#import stephane.display.graphes as graphes
import stephane.display.graphes as graphes

import garmin
import fitdecode
import tilemapbase

# Import the required library
from geopy.geocoders import Nominatim


def tmp_connect():
    tilemapbase.start_logging()
    tilemapbase.init(create=True)
    # Use open street map
    t = tilemapbase.tiles.build_OSM()
    return t

def project(Long,Lat):
    Long = np.asarray(Long)
    Lat = np.asarray(Lat)
    
    xtile = (Long + 180.0) / 360.0
    lat_rad = np.radians(Lat)
    
    ytile = (1.0 - np.log(np.tan(lat_rad) + (1 / np.cos(lat_rad))) / np.pi) / 2.0
    return (xtile, ytile)

def extent(BBox):
    #to be revised, weird projection system
    ext = tilemapbase.Extent.from_lonlat(BBox[0],BBox[1],BBox[2],BBox[3])
    ext = ext.to_aspect(1.0)
    return ext

def display_map(ext,t,title='',ax=None,width=600):
    plotter = tilemapbase.Plotter(ext, t, width=width)
    plotter.plot(ax, t)

    figs = {}
    figs.update(graphes.legende('Longitude','Latitude',title,cplot=True))

    return ax,figs

def map_traj(Long,Lat,save=False,scale=1.2,title=''):
    BBox = box_data(Long,Lat,scale=scale)
    print(BBox)
    X,Y = project(Long,Lat)
    ext = extent(BBox)

    t = tmp_connect()
    
    fig, ax = plt.subplots(figsize=(8, 8), dpi=200)
    ax,figs = display_map(ext,t,ax=ax,width=600)
    ax.plot(X,Y,'k-')

    figs.update(graphes.legende('Longitude','Latitude',title,cplot=True))

    if save:
        graphes.save_figs(figs,savedir=savefolder,frmt='png')

    return ax,t,figs


def box_data(Longs,Lats,scale=1.2,width=0.02,square=True):
    #return a bounding box scaled 
    X,Y = Longs,Lats

    Wlong = (np.max(X)-np.min(X))*scale
    Wlat = (np.max(Y)-np.min(Y))*scale
    Xc = (np.max(X)+np.min(X))/2
    Yc = (np.max(Y)+np.min(Y))/2

    if square:
        W = np.max([Wlong,Wlat])
        Wlong = W
        Wlat = W
    BBox = [Xc-Wlong,Xc+Wlong,Yc-Wlat,Yc+Wlat]
    #bound the BBox
    BBox[0] = np.max([BBox[0],-179.9])
    BBox[1] = np.min([BBox[1],179.9])
    BBox[2] = np.max([BBox[2],-89.9])
    BBox[3] = np.min([BBox[3],89.9])    
    return BBox

def boxes(name):
    d={}
    d['haha'] = [-68.8311,-68.7995,48.3389,48.3551]
    d['bic'] = [-68.8681,-68.7995,48.3181,48.3551]
    d['rimouski'] = [-68.5462,-68.5065,48.4374,48.4636]
    d['capelans'] = [-68.8625,-68.8330,48.3185,48.3363]
    return d[name]

def check_box(Long,Lat,BBox):
#    BBox = boxes(name)
    bLong = np.logical_and(Long>BBox[0],Long<BBox[1])
    bLat = np.logical_and(Lat>BBox[2],Lat<BBox[3])
    b= np.logical_and(bLong,bLat)
    print('Box dimension :'+str(b))
    return b

def get_city(Long,Lat):
# Initialize Nominatim API
    Longc,Latc = np.mean(Long),np.mean(Lat)
    geolocator = Nominatim(user_agent="MyApp")
    location = geolocator.reverse((Latc,Longc))
    return location

def display_city(Long,Lat):
    location = get_city(Long,Lat)
    print(location)
    if location is None:
        return ''
    parse = location.address.split(',')
    s = parse[2][1:]+','+parse[-1] #fucking arbitrary
    print(s)
    return s

def title_std(filename,Long,Lat):
    date = filename.split('/')[-1].split(' ')[0]
    hour = filename.split('/')[-1].split(' ')[1].split('.fit')[0]
    print(date,hour)

    s = display_city(Long,Lat)
    title = date+', '+hour.replace('.','_')+', '+s
    return title

    
    
