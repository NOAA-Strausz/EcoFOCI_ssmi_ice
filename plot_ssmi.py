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
from matplotlib.mlab import griddata
import cmocean
import sys

data_file=sys.argv[1]
#data_file='nt_20180402_f18_nrt_n.bin'
latfile='/home/dave/workstuff/github/EcoFOCI_ssmi_ice/psn25lats_v3.dat'
lonfile='/home/dave/workstuff/github/EcoFOCI_ssmi_ice/psn25lons_v3.dat'

def decode_datafile(filename):
    icefile = open(filename, 'rb')
    #remove the header
    icefile.seek(300)
    ice = np.fromfile(icefile,dtype=np.uint8)
    ice[ice >= 253] = 0
    ice = ice/2.5
    return ice;

def get_date(filename):
    icefile = open(filename, 'rb')
    header = icefile.read(300)
    #the date is located btw the following bytes
    date = date=header[219:229].decode() #gets date from header
    date = dt.datetime.strptime(date,"%m/%d/%Y")
    return date;

def decode_latlon(filename):
    latlon_file = open(filename, 'rb')
    output = np.fromfile(latlon_file,dtype='<i4')
    output = output/100000.0
    return output;

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
#land_50m = cfeature.NaturalEarthFeature('physical', 'land', '50m',
#                                        edgecolor='face',
#                                        facecolor=cfeature.COLORS['land'])

bath_50m = cfeature.NaturalEarthFeature('raster', 'OB_50M', '50m')
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

extent = [160, 222, 50, 75]
projection=ccrs.LambertConformal(central_longitude=200.0)
transformation=ccrs.PlateCarree()
fig, ax = make_map(projection=projection)

#ax.plot(df.longitude.values,df.latitude.values,'k.',markersize=.25,transform=transformation)
cm=ax.pcolormesh(xi,yi,zi,transform=transformation,cmap=cmocean.cm.ice, vmin=0, vmax=100)
plt.colorbar(cm)
ax.stock_img()
#ax.add_feature(land_50m)
#ax.add_feature(bath_50m)
ax.coastlines(resolution='50m')
ax.set_extent(extent)
plot_title="Ice Concentration: "+file_date.strftime("%Y-%m-%d")
t = fig.suptitle(plot_title)
filename=filename_prefix + '_ice_plot'
fig.savefig(filename)
