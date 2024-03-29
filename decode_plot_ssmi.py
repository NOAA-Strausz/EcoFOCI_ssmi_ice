#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  8 14:48:21 2019

@author: dave
"""

import numpy as np
import pandas as pd
import datetime as dt
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.io import shapereader
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from scipy import interpolate 
#from matplotlib.mlab import griddata
import cmocean
import sys
import argparse
import matplotlib.ticker as mticker
import cartopy.mpl.ticker as cticker
from datetime import datetime

def valid_date(s):
    try:
        return datetime.strptime(s, "%Y-%m-%d")
    except ValueError:
        msg = "Not a valid date: '{0}'.".format(s)
        raise argparse.ArgumentTypeError(msg)

parser = argparse.ArgumentParser(description='Decode binary SSMI satellite data')
parser.add_argument('infile', nargs='?',type=str,
                    help='full path to input file')
parser.add_argument('-d', '--date', nargs=1, type=valid_date,
                    help='pull file for particular date. do ont use with infile, date must be Y-m-d')
parser.add_argument('-m', '--mooring', nargs=1, 
                    help='add mooring location, ie ck1-ck4 or bs2-bs8 or all')
parser.add_argument('-ex', '--extents', nargs=1,
                    help='chooses extents of map, options are bering, chukchi, custom, and default')
parser.add_argument('-p', '--plot', action="store_true",
                    help='plot the ice concentrations on a map')
parser.add_argument('-c', '--csv', action="store_true",
                    help='make a CSV file of the data with lat,lon, and ice percent')


args=parser.parse_args()




#data_file=sys.argv[1]
#data_file='nt_20180402_f18_nrt_n.bin'
#these are the files that contain the lats and lons. Obtained from here:
#ftp://sidads.colorado.edu/pub/DATASETS/seaice/polar-stereo/tools/
latfile='/home/makushin/strausz/ecofoci_github/EcoFOCI_ssmi_ice/psn25lats_v3.dat'
lonfile='/home/makushin/strausz/ecofoci_github/EcoFOCI_ssmi_ice/psn25lons_v3.dat'
nrt_path='/home/akutan/strausz/ssmi_ice/data/nrt/'
boot_path='/home/akutan/strausz/ssmi_ice/data/bootstrap/'

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
    #output = int(output * 1000)/1000 #sets decimal place at 3 without rounding
    return output;

if args.infile:        
    data_file=args.infile
    data={'latitude':decode_latlon(latfile), 'longitude':decode_latlon(lonfile),
          'ice_conc':decode_datafile(data_file)}
    df=pd.DataFrame(data)
    
    file_date=get_date(data_file)
    filename_prefix=file_date.strftime("%Y_%m_%d")
elif args.date:
    #for now only works with nrt
    data_file=nrt_path + 'nt_' + args.date[0].strftime("%Y%m%d") + '_f18_nrt_n.bin'
    data={'latitude':decode_latlon(latfile), 'longitude':decode_latlon(lonfile),
          'ice_conc':decode_datafile(data_file)}
    df=pd.DataFrame(data)
    
    filename_prefix=args.date[0].strftime("%Y_%m_%d")
    file_date = args.date[0]
    
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
    fig, ax = plt.subplots(figsize=(10.5, 7),
                           subplot_kw=dict(projection=projection))
    if projection == ccrs.PlateCarree():
        gl = ax.gridlines(draw_labels=True)
        gl.xlabels_top = gl.ylabels_right = False
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
    return fig, ax

#download land mask
# 50m is a good balance between dataset size and land feature resolution

#if csv file is wanted
if args.csv:
    filename = filename_prefix + "_ice_conc.csv"
    df.to_csv(filename)
    
    
if args.plot:
    land_50m = cfeature.NaturalEarthFeature('physical', 'land', '50m',
                                            edgecolor='face',
                                            facecolor=cfeature.COLORS['land'])
    
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
    moorings = {'ck1': [70.838,163.125], 'ck2': [71.231,164.223], 'ck3': [71.828,166.070],
                'ck4': [71.038,160.514], 'bs2': [56.869,164.050], 'bs4': [57.895,168.878],
                'bs5': [59.911,171.731], 'bs8': [62.194,174.688]}
    
    
    
    #dictionary for various extents of map
    
    extent = {'chukchi': [180, 210, 64, 73], 'bering': [174,206,51,66], 
              'default': [165,220,50,75],'custom':[180,228,68.25,77.25]} #adjust custom to desired extents
    
    if args.extents:
        map_extents=extent[args.extents[0]]
    else:
        map_extents=extent['default']
        
    #projection=ccrs.LambertConformal(central_longitude=200.0)
    projection=ccrs.Mercator(central_longitude=180)
    #projection=ccrs.PlateCarree(central_longitude=200.0)
    transformation=ccrs.PlateCarree()
    fig, ax = make_map(projection=projection)
    
    cm=ax.pcolormesh(xi,yi,zi,transform=transformation,cmap=cmocean.cm.ice, vmin=0, vmax=100)
    plt.colorbar(cm)
    
    ax.add_feature(land_50m)
    ax.coastlines(resolution='50m')
    plot_title="Ice Concentration: "+file_date.strftime("%Y-%m-%d")
    
    if args.mooring:
        if args.mooring[0] == 'all':
            for key in moorings:
                lat = moorings[key][0]
                lon = 360-moorings[key][1]
                ax.plot(lon, lat, 'r.', transform=transformation)
                #plot_title = key.swapcase() + ' ' + plot_title
                label = key.swapcase()
                label_lat = lat
                label_lon = lon
                ax.text(label_lon, label_lat, label, horizontalalignment='right', 
                        verticalalignment='bottom', transform=transformation)
                
        else:    
            lat = moorings[args.mooring[0]][0]
            lon = 360-moorings[args.mooring[0]][1]
            ax.plot(lon, lat, 'r.', transform=transformation)
            plot_title = args.mooring[0].swapcase() + ' ' + plot_title
            label = args.mooring[0].swapcase()
            label_lat = lat
            label_lon = lon
            ax.text(label_lon, label_lat, label, horizontalalignment='right', 
                    verticalalignment='bottom', transform=transformation)
    
    
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
