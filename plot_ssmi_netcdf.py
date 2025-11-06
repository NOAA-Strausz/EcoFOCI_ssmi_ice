#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  8 14:48:21 2019
                    
@author: dave
"""

import numpy as np
import pandas as pd
import xarray as xr
import re
import datetime as dt
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.io import shapereader
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from scipy import interpolate 
from scipy.interpolate import griddata
import cmocean
import sys
import argparse
import matplotlib.ticker as mticker
import cartopy.mpl.ticker as cticker
import yaml
import pathlib

parser = argparse.ArgumentParser(description='Decode binary SSMI satellite data')
parser.add_argument('infile', metavar='infile', type=str,
                    help='full path to input file')
parser.add_argument('-m', '--mooring', nargs='+', 
                    help='add mooring location, individual moorings or [all] for all moorings \
                    must be after the infile if listing individual moorings [ie -m c1 c2 m2 m4]')
parser.add_argument('-ex', '--extents', nargs=1,
                    help='chooses extents of map, options are bering, chukchi, custom, and default')


args=parser.parse_args()

#this gets the path of where the executable file is located so that you can
#run from anywhere and don't have to tell it where the config file is
file_path=str(pathlib.Path(__file__).parent.resolve())
config_file=file_path + '/' + 'ice_config.yaml'
#get config settings from yaml file
#view yaml setup file for descriptions of these variables
with open(config_file, 'r') as file:
    config = yaml.safe_load(file)


#data_file=sys.argv[1]
#data_file='nt_20180402_f18_nrt_n.bin'
#these are the files that contain the lats and lons. Obtained from here:
#ftp://sidads.colorado.edu/pub/DATASETS/seaice/polar-stereo/tools/
latfile = config['latfile']
lonfile = config['lonfile']
boot_year = config['boot_year']

def decode_datafile(filename):
    ds = xr.open_dataset(filename)
    #search for the variable with "ICECON" in it
    #bootstrap files should only have one data variable so the below
    #code searches for the variable with 'ICECON' in the name
    #it changes over the years depending on which sattelite
    #the nrt files have three data variables, F16_ICECON, F17_ICECON, 
    #and F18_ICECON.  We will just use F18 for now. 
    variable=[var for var in ds.data_vars if 'ICECON' in var]
    variable=variable[-1:].pop() #gets the last variable in the list
    ice = ds[variable].values.flatten()
    # if re.search("NSIDC0081", filename):
        # variable=[var for var in ds.data_vars if 'ICECON' in var]
        # ice = ds[variable].values.flatten()
    # else:
        # ice = ds['F18_ICECON'].values.flatten()
    ice = ice*100
    # ~ ice[ice >= 253] = 0
    # ~ ice = ice/2.5
    # ~ ice = np.round(ice, 1)
    #below is old code for decoding the binary files
    #not needed since 2024 when bootstrap went to all netcdf
    # else:
        # prefix = filename.split('/')[-1:][0][:2] 
        # icefile = open(filename, 'rb')
    
        # if prefix == 'nt':
            # #remove the header
            # icefile.seek(300)
            # ice = np.fromfile(icefile,dtype=np.uint8)
            # ice[ice >= 253] = 0
            # ice = ice/2.5
        # elif prefix == 'bt':
            # ice = np.fromfile(icefile,dtype=np.uint16)
            # ice = ice/10.
            # ice[ice == 110] = 100 #110 is polar hole
            # ice[ice == 120] = np.nan #120 is land
        # else: 
            # ice=np.nan
        
        # icefile.close()
    
    return ice;

def get_date(filename):
    #gets date from filename
    #first remove path from filename if it is there
    filename = filename.split('/')[-1:][0]
    date = re.search("\\d{8}", filename).group()        
    date = dt.datetime.strptime(date,"%Y%m%d")
    return date;

def decode_latlon(filename):
    latlon_file = open(filename, 'rb')
    output = np.fromfile(latlon_file,dtype='<i4')
    output = output/100000.0
    #output = int(x * 1000)/1000 #sets decimal place at 3 without rounding
    return output;

if args.infile:        
    data_file=args.infile
    data={'latitude':decode_latlon(latfile), 'longitude':decode_latlon(lonfile),
          'ice_conc':decode_datafile(data_file)}
    df=pd.DataFrame(data)
    
    file_date=get_date(data_file)
    filename_prefix=file_date.strftime("%Y_%m_%d")

### set a range of lats and lons
# not advised as this messes with gridding assumptions later
#df.drop(df.loc[((df['latitude']<=45) | (df['latitude']>=75))].index, inplace=True)
#df.drop(df.loc[((df['longitude']<=-180) | (df['longitude']>=-150))].index, inplace=True)

### remove 0's by either dropping and making database smaller or by replacing with nans
# not advised for highly regional views as code may interpolate over these gaps oddly 
# . or zeros may be actual (polynas).  Global views may be ok and speed up though

#df.drop(df.loc[df['ice_conc']==0].index, inplace=True)
# or
df['ice_conc'][df['ice_conc']==0] = np.nan
#df['ice_conc'][df['ice_conc']<100] = np.nan
#df['ice_conc'][df['ice_conc']<0] = np.nan
def make_map(projection=ccrs.PlateCarree()):
    fig, ax = plt.subplots(figsize=(7, 7),
                           subplot_kw=dict(projection=projection))
    if projection == ccrs.PlateCarree():
        gl = ax.gridlines(draw_labels=True)
        gl.xlabels_top = gl.ylabels_right = False
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
    return fig, ax

#download land mask
# 50m is a good balance between dataset size and land feature resolution
land_50m = cfeature.NaturalEarthFeature('physical', 'land', '50m',
                                        edgecolor='face',
                                        facecolor=cfeature.COLORS['land'])
#for zoomed in, more detail
land_10m = cfeature.NaturalEarthFeature('physical', 'land', '10m',
                                        edgecolor='face',
                                        facecolor=cfeature.COLORS['land'])

#bath_50m = cfeature.NaturalEarthFeature('raster', 'OB_50M', '50m')
#                                        edgecolor='face',
#                                        facecolor=cfeature.COLORS['water'])


### Plot location of sample points for general BS region
#extent = [180, 210, 55, 68]
#projection=ccrs.LambertConformal(central_longitude=200.0)
#transformation=ccrs.PlateCarree()
#fig,ax = make_map(projection=projection)
#
#ax.plot(df.longitude,df.latitude,'k.',markersize=.25,transform=transformation)
#ax.add_feature(land_50m)
#ax.coastlines(resolution='50m')
#ax.set_extent(extent)

### Remapping the modes from the analysis
#-- Now let's grid your data.
# First we'll make a regular grid to interpolate onto. This is equivalent to
# your call to `mgrid`, but it's broken down a bit to make it easier to
# understand. 

#The number of columns and rows can be directly linked to the grid resolution
#360 cols would be 1deg resolution... 180 cols would be 2deg resolution
#the more columns, the slower the gridding process but the smoother the plot.  Too
#many columns will lead to oversampling so .25x.25 is probably the highest I would go

numcols, numrows = 360*4, 90*4

xi = np.linspace(df.longitude.min(), df.longitude.max(), numcols)
yi = np.linspace(df.latitude.min(), df.latitude.max(), numrows)
xi, yi = np.meshgrid(xi, yi)

#-- Interpolate at the points in xi, yi
# "griddata" expects "raw" numpy arrays, so we'll pass in
# data.x.values instead of just the pandas series data.x

#%%timeit

# regridding data with 0's removed data
x, y, z = df.longitude.values, df.latitude.values, df.ice_conc.values
zi = interpolate.griddata((x, y),z, (xi, yi), method='linear')

#adds point on map
#dictinary of mooring locations
moorings = {'c1': [70.838,163.125], 'c2': [71.231,164.223], 'c3': [71.828,166.070],
            'c4': [71.038,160.514], 'c5': [71.203,158.011], 'c12': [67.911,168.195],
            'm2': [56.869,164.050], 'm4': [57.895,168.878], 'm5': [59.911,171.731],
            'm8': [62.194,174.688], 'm14': [64.0,167.929], 'k1': [57.855,163.667],
            'k2': [58.262,163.495], 'k3': [58.785,163.278]}



#dictionary for various extents of map
#note that longitude is positive degrees west

extent = {'chukchi': [180, 210, 64, 73], 'bering': [174,206,51,66], 
          'default': [165,220,50,75],'custom':[185,200,56,61]} #adjust custom to desired extents

if args.extents:
    map_extents=extent[args.extents[0]]
else:
    map_extents=extent['default']
    
#projection=ccrs.LambertConformal(central_longitude=200.0)
projection=ccrs.Mercator(central_longitude=180)
#projection=ccrs.PlateCarree(central_longitude=200.0)
transformation=ccrs.PlateCarree()
fig, ax = make_map(projection=projection)

#ax.plot(df.longitude.values,df.latitude.values,'k.',markersize=.25,transform=transformation)
cm=ax.pcolormesh(xi,yi,zi,transform=transformation,cmap=cmocean.cm.ice, vmin=0, vmax=100)
plt.colorbar(cm)
#ax.stock_img()
ax.add_feature(land_50m)
#uncomment out for more detail, will be slower
#ax.add_feature(land_10m)
ax.coastlines(resolution='50m')
#uncomment out for more detail, will be slower
#ax.coastlines(resolution='10m')

plot_title="Ice Concentration: "+file_date.strftime("%Y-%m-%d")
#uncomment out below for custom plot title
# ~ plot_title="Bering Sea Maximum Ice Extent, Mar 19, 2024"
if args.mooring:
    if args.mooring[0] == 'all':
        selected_moorings = list(moorings.keys())
    else:
        selected_moorings = args.mooring
    for i in selected_moorings:
        lat = moorings[i][0]
        lon = 360-moorings[i][1]
        ax.plot(lon, lat, 'r.', transform=transformation)
        #plot_title = key.swapcase() + ' ' + plot_title
        label = i.swapcase()
        label_lat = lat
        label_lon = lon
        ax.text(label_lon, label_lat, label, horizontalalignment='right', 
                verticalalignment='bottom', transform=transformation)
                
#below area can be used to add custom annotations to the map:
date = file_date.strftime("%b %d")
ax.text(210, 52, date, horizontalalignment='right', 
                verticalalignment='bottom', transform=transformation,
                size='x-large', weight='bold')
            
#attempt at making gridlines with cartopy stuff, doesn't work that well
#gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
#                  linewidth=1, color='gray', alpha=0.5, linestyle='-')
#gl.xlabels_top = False
#gl.xlabels_bottom = True
#gl.ylabels_right = False
#gl.xlines = True
#gl.ylines = True
#gl.xformatter = LONGITUDE_FORMATTER
#gl.yformatter = LATITUDE_FORMATTER


        
#ax.set_xticks([178, -178, -174, -170, -166, -162], crs=ccrs.PlateCarree())
#ax.set_xticklabels([178, -178, -174, -170, -166, -162])
#ax.set_yticks([55, 57, 59, 61, 63, 65], crs=ccrs.PlateCarree())
#ax.set_yticklabels([55, 57, 59, 61, 63, 65])
#lon_formatter = cticker.LongitudeFormatter()
#lat_formatter = cticker.LatitudeFormatter()
#ax.xaxis.set_major_formatter(lon_formatter)
#ax.yaxis.set_major_formatter(lat_formatter)
#ax.grid(linewidth=1, color='black', alpha=0.5, linestyle='-')
        

#use following to make ticks and grid if needed
#ax.set_xticks([-174, -168, -162, -156, -150, -144, -138, -132], crs=ccrs.PlateCarree())
#ax.set_xticklabels([-174, -168, -162, -156, -150, -144, -138, -132])
#ax.set_yticks([69.0, 70.5, 72.0, 73.5, 75.0, 76.5], crs=ccrs.PlateCarree())
#ax.set_yticklabels([69.0, 70.5, 72.0, 73.5, 75.0, 76.5])
#lon_formatter = cticker.LongitudeFormatter()
#lat_formatter = cticker.LatitudeFormatter()
#ax.xaxis.set_major_formatter(lon_formatter)
#ax.yaxis.set_major_formatter(lat_formatter)
#ax.grid(linewidth=1, color='black', alpha=0.5, linestyle='-')
        
ax.set_extent(map_extents)
t = fig.suptitle(plot_title)
filename=filename_prefix + '_ice_plot'
fig.savefig(filename)

#gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
#                  linewidth=2, color='gray', alpha=0.5, linestyle='--')
#gl.xlabels_top = False
#gl.xlabels_bottom = False
#gl.ylabels_right = False
#gl.xlines = False
#gl.ylines = False
#gl.xformatter = LONGITUDE_FORMATTER
#gl.yformatter = LATITUDE_FORMATTER
