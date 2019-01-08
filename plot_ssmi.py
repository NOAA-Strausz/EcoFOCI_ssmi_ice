#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  8 14:48:21 2019

@author: dave
"""

import numpy as np
import pandas as pd
import datetime as dt

data_file='nt_20180402_f18_nrt_n.bin'
latfile='psn25lats_v3.dat'
lonfile='psn25lons_v3.dat'

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
