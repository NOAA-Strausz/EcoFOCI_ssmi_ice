#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 20 13:15:11 2019

@author: strausz
"""

import argparse
import glob
import os
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
import warnings

parser = argparse.ArgumentParser(description='Get ice concentration around a point')
parser.add_argument('-latlon', '--latlon', nargs=2, 
                    help='latitude and longitude of desired point, W lon must be negative', type=float)
parser.add_argument('-y', '--years', nargs='+', help='year range ie "2015 2019"; can be a single year. Select "latest" to just print out latest data', 
                    type=int)
parser.add_argument('-l', '--latest', help='output latest value to stdout', action="store_true")
parser.add_argument('-L', '--linear', help='change output to one file with full time series in one column', action="store_true")
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

#suppress some warnings
warnings.filterwarnings("ignore")
#this gets the path of where the executable file is located so that you can
#run from anywhere and don't have to tell it where the config file is
file_path=str(pathlib.Path(__file__).parent.resolve())
file_path=file_path + '/'
config_file=file_path + 'meaniceinbox_config.yaml'
#get config settings from yaml file
#view yaml setup file for descriptions of these variables
with open(config_file, 'r') as file:
    config = yaml.safe_load(file)

latlon_nc = file_path + config['latlon_nc']
latfile = config['latfile']
lonfile = config['lonfile']
bootstrap = config['bootstrap']
nrt = config['nrt']
boot_year = config['boot_year']
amsr = config['amsr']
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
    date = pd.Timestamp(ds.time.values[0])
    return ice, date;

def decode_latlon(filename):
    ds = xr.open_dataset(filename)
    latitude = ds.latitude.values.flatten()
    longitude = ds.longitude.values.flatten()
    output = {'latitude': latitude, 'longitude': longitude}
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

def get_latest_file(directory):
    files = os.listdir(directory)
    full_paths = [os.path.join(directory, f) for f in files]
    latest_file = max(full_paths, key=os.path.getmtime)
    return latest_file
    
if args.radius:
    nlat, slat, wlon, elon = find_box(inlat, inlon, 2*args.distance, nm=args.nm)
    dist_type = 'Radius'
else:
    nlat, slat, wlon, elon = find_box(inlat, inlon, args.distance, nm=args.nm)
    dist_type = 'Box'
years = args.years
#put desired years in list
years = list(range(args.years[0],args.years[1]+1))
files = []

if args.latest:
    latest_file = get_latest_file(amsr)
    print("Newest file is " + latest_file)
    files.append(latest_file)
else:
    for i in years:
        year = str(i)
        if i <= boot_year:
            path = bootstrap + year + '/'
            file_list = glob.glob(path + '*.nc')
            file_list = sorted(file_list)
            files.extend(file_list)
        else:
            #changed in jan 2023 to only work with new netcdf nrt files
            path = amsr
            file_list = glob.glob(path + '*' + year + ('[0-9]' *4) + '*.nc')
            file_list = sorted(file_list)
            files.extend(file_list)

output_date = []
output_ice = []

#get latitude and longitude arrays from netcdf file
latlon = decode_latlon(latlon_nc)
for i in files:
    try:
        #print('decoding filename: ' + i)
        data_ice = latlon
        data_ice['ice_conc'], date = decode_datafile(i)
        df_ice=pd.DataFrame(data_ice)
        df_ice_chopped = df_ice[(df_ice.latitude <= nlat) & (df_ice.latitude >= slat) & 
                                (df_ice.longitude >= wlon) & (df_ice.longitude <= elon)]
        #multiply ice concentration by 100 to get it to percentage
        #xarray should have applied the mask_and_scale=True by default
        #when processing the nrt or bootstrap netcdf files to return a 
        #concentration that is from 0-1 where 1 is 100% ice cover
        
        df_ice_chopped.loc[:, 'ice_conc'] *= 100
        if args.radius:
            fn = lambda x: haversine((inlat, inlon),(x.latitude,x.longitude))
            distance = df_ice_chopped.apply(fn, axis=1)
            df_ice_chopped = df_ice_chopped.assign(dist=distance.values)
            df_ice_chopped = df_ice_chopped.loc[df_ice_chopped.dist < args.distance]        
        #date_string = date.strftime("%Y,%j")
        #round entire ice_chopped data frame to one decimal place
        df_ice_chopped=df_ice_chopped.round(1)
        ice = df_ice_chopped.ice_conc.mean().round(1)
        #below line worked until nans started showing up when only one point was available
        #ice = df_ice_chopped.ice_conc.mean().round(decimals=1)
        #print(date_string+','+str(ice))
        if args.verbose:
            print("Working on File: " + i)
            print("For " + date.strftime("%Y-%m-%d") + " ice concentration was " 
                  + str(ice) + "%")
        output_date.append(date)
        output_ice.append(ice)
    except:
        pass

data = {'date':output_date, 'ice_concentration': output_ice}
df = pd.DataFrame(data)
#rount to one decimal place for neatness
#not sure why it needs this again as it should be rounded already
#but nrt data was not getting rounded if multiple years selected
df['ice_concentration'] = df.ice_concentration.round(1)
df.set_index(['date'], inplace=True)
if args.linear or args.latest:
    df_out=df
else:
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
    


if args.latest:
    filename = ("meaniceinbox" + mooring + pointname + str(inlat) + lat_suffix + "_" + 
            str(abs(inlon)) + lon_suffix + "_" + str(args.distance) + units + "_" + dist_type + "_"
            + "latest_" + df_out.index[0].strftime("%Y-%m-%d") + ".txt")
    df_out.to_csv(filename, header=False)
elif args.linear:
    filename = ("meaniceinbox" + mooring + pointname + str(inlat) + lat_suffix + "_" + 
            str(abs(inlon)) + lon_suffix + "_" + str(args.distance) + units + "_" + dist_type + "_"
            + str(years[0]) + "-" + str(years[1]) + "_linear.csv")
    df_out.to_csv(filename)
else:
    if len(years) > 1:
        filename = ("meaniceinbox_" + mooring + pointname + str(inlat) + lat_suffix + "_" + 
                str(abs(inlon)) + lon_suffix + "_" + str(args.distance) + units + 
                "_" + dist_type + "_" + str(years[0]) + "-" + str(years[1]) + ".csv")
    else:
        filename = ("meaniceinbox_" + mooring + pointname + str(inlat) + lat_suffix + "_" + 
                str(abs(inlon)) + lon_suffix + "_" + str(args.distance) + units + 
                "_" + dist_type + "_" + str(years[0]) + ".csv")
    
        
    df_out.to_csv(filename)
    
