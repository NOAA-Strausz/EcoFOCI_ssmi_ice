#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 31 12:38:35 2019

@author: strausz
"""

import numpy as np
import datetime as dt
import argparse
#from scipy import interpolate 
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.io import shapereader
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.pyplot as plt
import cmocean
from matplotlib.image import imread
import cartopy.io.shapereader as shpreader


latfile='/home/makushin/strausz/ecofoci_github/EcoFOCI_ssmi_ice/psn25lats_v3.dat'
lonfile='/home/makushin/strausz/ecofoci_github/EcoFOCI_ssmi_ice/psn25lons_v3.dat'
#locations of ice files
bootstrap = '/home/akutan/strausz/ssmi_ice/data/bootstrap/'
nrt = '/home/akutan/strausz/ssmi_ice/data/nrt/'
background = '/home/makushin/strausz/map_data/OB_LR/OB_LR.tif'
shpfilename = '/home/makushin/strausz/map_data/ne_10m_bathymetry_J_1000/ne_10m_bathymetry_J_1000.shp'


parser = argparse.ArgumentParser(description='Plot ice on map of Alaska daily')
parser.add_argument('infile', metavar='infile', type=str, help='full path to input file')
args=parser.parse_args()

def decode_datafile(filename):
    #determine if it's nrt or bootstrap from filename prefix
    #note that we remove path first if it exists
    prefix = filename.split('/')[-1:][0][:2] 
    icefile = open(filename, 'rb')
    
    if prefix == 'nt':
        #remove the header
        icefile.seek(300)
        ice = np.fromfile(icefile,dtype=np.uint8)
        ice[ice >= 253] = 0
        ice = ice/2.5
    elif prefix == 'bt':
        ice = np.fromfile(icefile,dtype=np.uint16)
        ice = ice/10.
        ice[ice == 110] = 100 #110 is polar hole
        ice[ice == 120] = np.nan #120 is land
    else: 
        ice=np.nan
    
    icefile.close()
    
    return ice;

def get_date(filename):
    #gets date from filename
    #first remove path from filename if it is there
    filename = filename.split('/')[-1:][0]
    date = filename[3:11]
    date = dt.datetime.strptime(date,"%Y%m%d")
    return date;

def decode_latlon(filename):
    latlon_file = open(filename, 'rb')
    output = np.fromfile(latlon_file,dtype='<i4')
    output = output/100000.0
    #output = int(output * 1000)/1000 #sets decimal place at 3 without rounding
    latlon_file.close()
    return output;

def make_map(projection=ccrs.PlateCarree()):
    fig, ax = plt.subplots(figsize=(10.5, 7),
                           subplot_kw=dict(projection=projection))
    #if projection == ccrs.PlateCarree():
    #gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True)
#    gl.xlabels_top = gl.ylabels_right = False
#    gl.xformatter = LONGITUDE_FORMATTER
#    gl.yformatter = LATITUDE_FORMATTER
    return fig, ax, gl

if args.infile:        
    ice = decode_datafile(args.infile).reshape((448,304))
    ice[ice==0]=np.nan
    lat = decode_latlon(latfile).reshape((448,304))
    lon = decode_latlon(lonfile).reshape((448,304))
    
    x, y, z = lon.flatten(), lat.flatten(), ice.flatten()
    xsize = 448 * 4
    ysize = 448 * 4
    
    
    land_50m = cfeature.NaturalEarthFeature('physical', 'land', '50m',
                                            edgecolor='face',
                                            facecolor=cfeature.COLORS['land'])
    map_extents=[165,220,45,85]
    
    projection=ccrs.Orthographic(central_longitude=-170, central_latitude=90)
    transformation=ccrs.PlateCarree()
    fig, ax = plt.subplots(figsize=(12, 8),
                           subplot_kw=dict(projection=projection))
    cm=ax.pcolormesh(lon,lat,ice,transform=transformation,cmap=cmocean.cm.ice, vmin=0, vmax=100)
    #plt.colorbar(cm)
    
    #ax.add_feature(land_50m)
    ax.imshow(imread(background), origin='upper', transform=transformation,
              extent=[-180, 180, -90, 90])
    ax.coastlines(resolution='50m')
    gl = ax.gridlines(xlocs=range(0,390,10), ylocs=range(0,90,10))
    gl.n_steps=500
    ax.set_extent(map_extents)

    #ax.background_img(name='BM', resolution='high')
    file_date=get_date(args.infile)
    plot_title="Ice Concentration: "+file_date.strftime("%Y-%m-%d")
    t = fig.suptitle(plot_title)
    filename_prefix=file_date.strftime("%Y_%m_%d")
    filename=filename_prefix + '_ice_plot'
    fig.savefig(filename)
    