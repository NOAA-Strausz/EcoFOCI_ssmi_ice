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
latfile='psn25lats_v3.dat'
lonfile='psn25lons_v3.dat'

def decode_datafile(filename):
    icefile = open(filename, 'rb')
    ice = np.fromfile(icefile,dtype=np.uint8)
    ice[ice >= 253] = 0
    ice = ice/2.5
    return ice;

def get_date(filename):
    date = filename[3:10]
    date = dt.datetime.strptime(date,"%Y%m%d")
    
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

