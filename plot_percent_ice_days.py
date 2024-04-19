#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  11 2021
                    
@author: dave
"""

import argparse
import pathlib
import glob
import numpy as np
import yaml
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from scipy import interpolate
from scipy.interpolate import griddata
import cmocean

parser = argparse.ArgumentParser(description='Find amount of days in a year over certain ice %')
parser.add_argument('-y', '--year', nargs=1, help='Year that you want to process', 
                    type=int, required=True)
parser.add_argument('-m', '--mooring', nargs=1, 
                    help='add mooring location, ie ck1-ck4 or bs2-bs8 or all')
parser.add_argument('-ex', '--extents', nargs=1,
                    help='chooses extents of map, options are bering, chukchi, custom, and default')


args=parser.parse_args()

#this gets the path of where the executable file is located so that you can
#run from anywhere and don't have to tell it where the config file is
file_path=str(pathlib.Path(__file__).parent.resolve())
config_file=file_path + '/' + 'meaniceinbox_config.yaml'
#get config settings from yaml file
#view yaml setup file for descriptions of these variables
with open(config_file, 'r') as file:
    config = yaml.safe_load(file)

latfile = config['latfile']
lonfile = config['lonfile']
bootstrap = config['bootstrap']
nrt = config['nrt']
boot_year = config['boot_year']
file.close()

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
        ice[ice == 110] = 0 #110 is land
        ice[ice == 120] = 100 #120 is polar hole
    else: 
        ice=np.nan
    
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
    #output = int(x * 1000)/1000 #sets decimal place at 3 without rounding
    return output;

def make_map(projection=ccrs.PlateCarree()):
    fig, ax = plt.subplots(figsize=(10.5, 7),
                           subplot_kw=dict(projection=projection))
    if projection == ccrs.PlateCarree():
        gl = ax.gridlines(draw_labels=True)
        gl.xlabels_top = gl.ylabels_right = False
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
    return fig, ax



#generate list of paths for entire year
if args.year:  
    files=[]
    int_year = args.year[0]   
    year = str(int_year)
    if int_year <= boot_year:
        path = bootstrap + year + '/'
        files = files + glob.glob(path + '*.bin')
    else:
        path = nrt
        files = files + glob.glob(path + '*' + year + '*.bin')
    #make array with 0's that is length of binary data files
    #decoded binary files have 136192 values
    days_percent = np.zeros((136192))
    for i in files:
        print("Working on file: ", i)
        ice_conc = decode_datafile(i)
        ice_conc_orig = ice_conc
        ice_conc[ice_conc<15]=0
        ice_conc[ice_conc>=15]=1
        days_percent=days_percent + ice_conc
    data_ice={'latitude':decode_latlon(latfile), 'longitude':decode_latlon(lonfile),
          'days_above_percent':days_percent}
    df=pd.DataFrame(data_ice)
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
    x, y, z = df.longitude.values, df.latitude.values, df.days_above_percent.values
    zi = interpolate.griddata((x, y),z, (xi, yi), method='linear')  

    extent = {'chukchi': [180, 210, 64, 73], 'bering': [174,206,51,66], 
              'default': [165,220,50,75],'custom':[180,228,68.25,77.25]} #adjust custom to desired extents
    
    if args.extents:
        map_extents=extent[args.extents[0]]
    else:
        map_extents=extent['default']
            
        
    projection=ccrs.Mercator(central_longitude=180)
    transformation=ccrs.PlateCarree()
    fig, ax = make_map(projection=projection)
    
    land_50m = cfeature.NaturalEarthFeature('physical', 'land', '50m',
                                        edgecolor='face',
                                        facecolor=cfeature.COLORS['land'])
    
    cm=ax.pcolormesh(xi,yi,zi,transform=transformation,cmap=cmocean.cm.ice, vmin=0, vmax=366)
    plt.colorbar(cm)
    ax.add_feature(land_50m)
    ax.coastlines(resolution='50m')
    plot_title="Days above 15% ice in year: "+year
    
    ax.set_extent(map_extents)
    t = fig.suptitle(plot_title)
    
