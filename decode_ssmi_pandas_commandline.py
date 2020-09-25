#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  3 10:57:17 2019

@author: dave
"""

import numpy as np
import pandas as pd
import datetime as dt

import argparse

parser = argparse.ArgumentParser(description="SSMI Ice Binary File")
parser.add_argument(
    "DataPath",
    metavar="DataPath",
    type=str,
    help="full path to binary ssmi ice conc file",
)

args = parser.parse_args()

filename = args.DataPath
latfile = "psn25lats_v3.dat"
lonfile = "psn25lons_v3.dat"
areafile = "psn25area_v3.dat"


def decode_datafile(filename, datasource="NRT"):
    icefile = open(filename, "rb")

    if datasource is "NRT":
        # remove the header
        icefile.seek(300)
        ice = np.fromfile(icefile, dtype=np.uint8)
        ice[ice >= 253] = 0
        ice = ice / 2.5
    elif datasource is "Bootstrap":
        # no header
        ice = np.fromfile(icefile, dtype=np.uint16)
        ice = ice / 10
        # following are ingrained data masks... they do not coincide with the documentation
        # https://nsidc.org/data/nsidc-0079/versions/3
        ice[ice == 110] = 100  # 110 is polar hole
        ice[ice == 120] = 0  # 120 is land
    else:
        ice = np.nan

    return ice


def get_date(filename):
    icefile = open(filename, "rb")
    header = icefile.read(300)
    # the date is located btw the following bytes
    date = date = header[219:229].decode()  # gets date from header
    date = dt.datetime.strptime(date, "%m/%d/%Y")
    return date


def decode_latlon(filename):
    latlon_file = open(filename, "rb")
    output = np.fromfile(latlon_file, dtype="<i4")
    output = output / 100000.0
    return output


def decode_area(filename):
    # decode the areafile
    #area_file = open(areafile, "rb")
    area = np.fromfile(areafile, dtype="<i4")
    area = area / 1000.0
    return area


data = {
    "latitude": decode_latlon(latfile),
    "longitude": decode_latlon(lonfile),
    "grid_area": decode_area(areafile),
    "ice_conc": decode_datafile(filename, "NRT"),
}
df = pd.DataFrame(data)

df.to_csv(filename.replace(".bin", ".csv"), index=False)
