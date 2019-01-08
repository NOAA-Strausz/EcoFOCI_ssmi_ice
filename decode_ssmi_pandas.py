#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  3 10:57:17 2019

@author: dave
"""

import numpy as np
import pandas as pd
import datetime as dt

filename='nt_20180402_f18_nrt_n.bin'
latfile='psn25lats_v3.dat'
lonfile='psn25lons_v3.dat'

#decode the datafile
icefile = open(filename, 'rb')
header = icefile.read(300) #reads first 300 bytes which are the header
#doy = doy=header[215:218].decode()
date = date=header[219:229].decode() #gets date from header
date = dt.datetime.strptime(date,"%m/%d/%Y")
file_prefix=date.strftime("%Y_%j")
ice = np.fromfile(icefile,dtype=np.uint8)
ice[ice >= 253] = 0
ice = ice/10

#decode the latfile
lats_file = open(latfile, 'rb')
lats = np.fromfile(latfile,dtype='<i4')
lats = lats/100000.0

#decode the lonfile
lons_file = open(lonfile, 'rb')
lons = np.fromfile(lonfile,dtype='<i4')
lons= lons/100000.0

data={'latitude':lats, 'longitude':lons, 'ice_conc':ice}
df=pd.DataFrame(data)

df.to_csv(file_prefix+'.csv', index=False)