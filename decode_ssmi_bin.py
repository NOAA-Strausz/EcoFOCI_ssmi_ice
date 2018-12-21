#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 25 15:48:39 2018

@author: strausz
"""

#following procedure laid out on this website:
#https://geoinformaticstutorial.blogspot.com/2014/02/reading-binary-data-nsidc-sea-ice.html

import gdal
import struct
import numpy as np
import sys
import matplotlib.pyplot as plt


filename = sys.argv[1]

print "Processing file: "
print filename


#height and width of files
#see https://nsidc.org/data/polar-stereo/ps_grids.html
width = 304
height = 448

#for this code, inspiration found at https://stevendkay.wordpress.com/category/python/
icefile = open(filename, 'rb')
header = icefile.read(300)
ice = np.fromfile(icefile,dtype=np.uint8)
#ice = ice.reshape(height,width)
ice = ice.reshape(136192,1)
np.savetxt("test.csv", ice)
#ice = ice/10
ice[ice >= 253] = 0
ice = ice/10
pixels_with_ice = np.count_nonzero(ice)



icefile.close()

latfile = open('psn25lats_v3.dat', 'rb')
lats = np.fromfile(latfile,dtype='<i4')
lats = lats/100000.0
lats = lats.reshape(136192,1)
#lats = lats.reshape(height,width)

latfile.close()

lonfile = open('psn25lons_v3.dat', 'rb')
lons = np.fromfile(lonfile,dtype='<i4')
lons = lons/100000.0
lons = lons.reshape(136192,1)
lonfile.close()

ice_conc=np.hstack((lats,lons,ice))
np.savetxt('test_csv.csv',ice_conc, delimiter=',', fmt='%.3f')
#np.savetxt('test_lats.csv', lats, delimiter=',')

#ice = np.ma.masked_greater(ice,100)
#plt.imshow(ice)


## unpack binary data into a flat tuple z
#s="%dB" % (int(width*height),)
#z=struct.unpack_from(s, contents, offset = 300)
#
#nsidc = np.array(z).reshape((448,304))
#nsidc = np.rot90(nsidc, 1)
#
#
##write data to geotiff
#driver = gdal.GetDriverByName("GTiff")
#outraster = driver.Create('test4.tif', height, width, 1, gdal.GDT_Int16 )
#raster = np.zeros((width, height), np.float) 
#outraster.GetRasterBand(1).WriteArray( raster )
#
#outband = outraster.GetRasterBand(1)
#
##Write to file     
#outband.WriteArray(nsidc)
#outband.FlushCache()
#
##Clear arrays and close files
#outband = None
#iceraster = None
#outraster = None
#outarray = None