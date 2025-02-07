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
import sys
from haversine import haversine
import yaml
import pathlib
import xarray as xr
import re

parser = argparse.ArgumentParser(description='Get ice concentration around a point')
parser.add_argument('-latlon', '--latlon', nargs=2, 
                    help='latitude and longitude of desired point, W lon must be negative', type=float)
parser.add_argument('-y', '--years', nargs=2, help='year range ie "2015 2019"; must have at least 2 years, can be the same year', 
                    type=int, required=True)
parser.add_argument('-a', '--name', help='optional name of point', type=str)
parser.add_argument('-d', '--distance', help='size of box around point in km', 
                    type=float, required=True)
parser.add_argument('-n', '--nm', help='use nautical miles instead of km',
                    action="store_true")
parser.add_argument('-r', '--radius', help='use distance as radius around point instead of box',
                    action="store_true")
parser.add_argument('-m', '--mooring', help='Mooring name, choose from ck1-9, or bs2-14')
parser.add_argument('-v', '--verbose', help='Print some details while processing files',
                    action="store_true")
                                        
args=parser.parse_args()

#this gets the path of where the executable file is located so that you can
#run from anywhere and don't have to tell it where the config file is
file_path=str(pathlib.Path(__file__).parent.resolve())
config_file=file_path + '/' + 'meaniceinbox_config_pre_netcdf.yaml'
#get config settings from yaml file
#view yaml setup file for descriptions of these variables
with open(config_file, 'r') as file:
    config = yaml.safe_load(file)

latfile = config['latfile']
lonfile = config['lonfile']
bootstrap = config['bootstrap']
nrt = config['nrt']
nrt_ver2 = config['nrt_ver2']
boot_year = config['boot_year']
file.close()

#mooring locaitons taken from 'https://www.pmel.noaa.gov/foci/foci_moorings/mooring_info/mooring_location_info.html'
moorings = {'bs2':[56.869,-164.050], 'bs4':[57.895,-168.878], 'bs5':[59.911,-171.73],
         'bs8':[62.194,-174.688], 'bs14':[64.00,-167.933], 'ck1':[70.838,-163.125], 'ck2':[71.231,-164.223],
         'ck3':[71.828,-166.070], 'ck4':[71.038,-160.514], 'ck5':[71.203,-158.011],
         'ck9':[72.464,-156.548], 'ck10':[70.211,-167.787],
         'ck11':[70.013,-166.855], 'ck12':[67.911,-168.195]}

if args.mooring:
    if args.mooring in moorings:
        inlat = moorings[args.mooring][0]
        inlon = moorings[args.mooring][1]
    else:
        sys.exit("Mooring not listed")
    mooring = args.mooring + "_"
else:
    inlat = args.latlon[0]
    inlon = args.latlon[1]
    mooring = ''

if args.nm:
    units = 'nm'
else:
    units = 'km'
    
if args.name:
    pointname = args.name + "_"
else:
    pointname = ''

def decode_datafile(filename):
    #determine if it's nrt or bootstrap from filename prefix
    #note that we remove path first if it exists
    if filename[-2:] == 'nc': #looks to see if it a netcdf file
        ds = xr.open_dataset(filename)
        ice = ds['F18_ICECON'].values.flatten()
        ice = ice*250
        ice[ice >= 253] = 0
        ice = ice/2.5
        ice = np.round(ice, 1)
    else:
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
    #filename = filename.split('/')[-1:][0]
    #date = filename[3:11]
    #use regex to get date from filename
    date=re.search(r'\d{8}', filename).group(0)
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

    
if args.radius:
    nlat, slat, wlon, elon = find_box(inlat, inlon, 2*args.distance, nm=args.nm)
    dist_type = 'Radius'
else:
    nlat, slat, wlon, elon = find_box(inlat, inlon, args.distance, nm=args.nm)
    dist_type = 'Box'
#put desired years in list
years = list(range(args.years[0],args.years[1]+1))
files = []

for i in years:
    year = str(i)
    if i <= boot_year:
        path = bootstrap + year + '/'
        files = files + glob.glob(path + '*.bin')
    else:
        #changed in jan 2023 to only work with new netcdf nrt files
        path = nrt_ver2
        files = files + glob.glob(path + '*' + year + ('[0-9]' *4) + '*.nc')
        files = sorted(files)

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
    if args.radius:
        fn = lambda x: haversine((inlat, inlon),(x.latitude,x.longitude))
        distance = df_ice_chopped.apply(fn, axis=1)
        df_ice_chopped = df_ice_chopped.assign(dist=distance.values)
        df_ice_chopped = df_ice_chopped.loc[df_ice_chopped.dist < args.distance]        
    #date_string = date.strftime("%Y,%j")
    #round entire ice_chopped data frame to one decimal place
    df_ice_chopped=df_ice_chopped.round(1)
    ice = df_ice_chopped.ice_conc.mean()
    #below line worked until nans started showing up when only one point was available
    #ice = df_ice_chopped.ice_conc.mean().round(decimals=1)
    #print(date_string+','+str(ice))
    if args.verbose:
        print("Working on File: " + i)
        print("For " + date.strftime("%Y-%m-%d") + " ice concentration was " 
              + str(ice) + "%")
    output_date.append(date)
    output_ice.append(ice)

data = {'date':output_date, 'ice_concentration': output_ice}
df = pd.DataFrame(data)
#rount to one decimal place for neatness
df['ice_concentration'] = df.ice_concentration.round()
df.set_index(['date'], inplace=True)
years_grouped = df.groupby(df.index.year)

#make empty dataframe to keep structure of 1..366 for doy
df_out=pd.DataFrame(index=range(1,366))
for name, group in years_grouped:
    year = str(name)
    group.index=group.index.dayofyear
    group.index.rename("day_of_year", inplace=True)
    group.rename(columns={"ice_concentration" : year}, inplace=True)

    df_out = pd.concat([df_out,group], axis=1)
    
    
#df_out = pd.DataFrame.from_dict(output, orient='index').transpose()
#df_out['DOY']=df_out.DOY.astype(int)
#get longiude suffix, assum lat is north
lat_suffix = 'N'

if inlon < 0:
    lon_suffix = 'W'
else:
    lon_suffix = 'E'
    

filename = ("meaniceinbox_" + mooring + pointname + str(inlat) + lat_suffix + "_" + 
            str(abs(inlon)) + lon_suffix + "_" + str(args.distance) + units + 
            "_" + dist_type + "_" + str(args.years[0]) + "-" + str(args.years[1]) + ".csv")
df_out.to_csv(filename)
