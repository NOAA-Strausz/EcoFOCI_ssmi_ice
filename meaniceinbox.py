#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 20 13:15:11 2019

@author: strausz
"""

import argparse
import glob
import numpy as np
import datetime as dt
import pandas as pd
import math



parser = argparse.ArgumentParser(description='Get ice concentration around a point')
parser.add_argument('-latlon', '--latlon', nargs=2, 
                    help='latitude and longitude of desired point', type=float)
parser.add_argument('-y', '--years', nargs=2, help='year range ie "2015 2019"', 
                    type=int)
parser.add_argument('-d', '--distance', nargs=1, help='size of box around point', 
                    type=float)
parser.add_argument('-n', '--nm', help='use nautical miles instead of km',
                    action="store_true")
parser.add_argument('-r', '--radius', help='use circle of radius 'r' instead of box',
                    action="store_true")
args=parser.parse_args()

latfile='/home/makushin/strausz/ecofoci_github/EcoFOCI_ssmi_ice/psn25lats_v3.dat'
lonfile='/home/makushin/strausz/ecofoci_github/EcoFOCI_ssmi_ice/psn25lons_v3.dat'
#locations of ice files
bootstrap = '/home/akutan/strausz/ssmi_ice/data/bootstrap/'
nrt = '/home/akutan/strausz/ssmi_ice/data/nrt/'
#latest available bootstrap year
boot_year = 2018

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

def find_box(lat1, lon1, dist, nm):
    
    #formula pulled from this website:
    #http://www.movable-type.co.uk/scripts/latlong.html
    
    if nm:
        r=3440
    else:
        r=6371
    
    dist = dist/2
    
    wlon = math.radians(lon1) + math.atan2(math.sin(math.radians(270)) *
                        math.sin(dist/r) * math.cos(math.radians(lat1)),
                        math.cos(dist/r) - math.sin(math.radians(lat1)) *
                        math.sin(math.radians(lat1)))
    elon = math.radians(lon1) + math.atan2(math.sin(math.radians(90)) *
                        math.sin(dist/r) * math.cos(math.radians(lat1)),
                        math.cos(dist/r) - math.sin(math.radians(lat1)) *
                        math.sin(math.radians(lat1)))
    nlat = math.asin(math.sin(math.radians(lat1)) * math.cos(dist/r) + 
                     math.cos(math.radians(lat1)) * math.sin(dist/r) * 
                     math.cos(math.radians(0)))
    slat = math.asin(math.sin(math.radians(lat1)) * math.cos(dist/r) + 
                     math.cos(math.radians(lat1)) * math.sin(dist/r) * 
                     math.cos(math.radians(180)))
    
    wlon = round(math.degrees(wlon), 4)
    elon = round(math.degrees(elon), 4)
    nlat = round(math.degrees(nlat), 4)
    slat = round(math.degrees(slat), 4)
    
    return([nlat,slat,wlon,elon])

#put desired years in list
    
nlat, slat, wlon, elon = find_box(args.latlon[0], 
                                  args.latlon[1], args.distance[0], nm=args.nm)
years = list(range(args.years[0],args.years[1]+1))
files = []

for i in years:
    year = str(i)
    if i <= boot_year:
        path = bootstrap + year + '/'
        files = files + glob.glob(path + '*.bin')
    else:
        path = nrt
        files = files + glob.glob(path + '*' + year + '*.bin')

output_date = []
output_ice = []
        
for i in files:
    #print('decoding filename: ' + i)
    data_ice={'latitude':decode_latlon(latfile), 'longitude':decode_latlon(lonfile),
          'ice_conc':decode_datafile(i)}
    df_ice=pd.DataFrame(data_ice)
    df_ice_chopped = df_ice[(df_ice.latitude <= nlat) & (df_ice.latitude >= slat) & 
                            (df_ice.longitude >= wlon) & (df_ice.longitude <= elon)]
    date = get_date(i)
    #date_string = date.strftime("%Y,%j")
    ice = df_ice_chopped.ice_conc.mean().round(decimals=1)
    #print(date_string+','+str(ice))
    output_date.append(date)
    output_ice.append(ice)

data = {'date':output_date, 'ice_concentration': output_ice}
df = pd.DataFrame(data)
df.set_index(['date'], inplace=True)
years_grouped = df.groupby(df.index.year)

dummy_list = []
output={'doy':range(1,367)}
for name, group in years_grouped:
    year = str(name)
    if name <= 1988:
        group = group.resample('d').mean()
        if name == 1978:
            dummy_list = [np.nan]*304
        elif name in [1979, 1982, 1984, 1985, 1987]:
            dummy_list = [np.nan]
        elif name == 1988:
            dummy_list = [np.nan]*13
        else:
            dummy_list = []
        output[year] = dummy_list + list(group.ice_concentration)
    else:
        output[year] = list(group.ice_concentration)
    
    
df_out = pd.DataFrame.from_dict(output, orient='index').transpose()
#get longiude suffix, assum lat is north
lat_suffix = 'N'
if args.years[1] < 0:
    lon_suffix = 'W'
else:
    lon_suffix = 'E'
filename = ("meaniceinbox_" + str(args.latlon[0]) + lat_suffix + "_" + 
            str(abs(args.latlon[1])) + lon_suffix + "_" + str(args.years[0]) + 
            "-" + str(args.years[1]) + ".csv")
df_out.to_csv(filename, index=False)