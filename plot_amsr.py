#!/usr/bin/env python3
#using this dataset:  https://nsidc.org/data/nsidc-0803/versions/2


import argparse
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from pyproj import Proj
import cartopy.feature as cfeature
import matplotlib.colors as colors
import pandas as pd
import matplotlib.ticker as plticker
import cmocean
import re
from datetime import datetime


parser = argparse.ArgumentParser(description='Plot NSIDC-0803, AMSR2 Daily Polar Gridded Sea Ice Concentrations, Version 2 Data Files ')
parser.add_argument('infile', metavar='infile', type=str,
                    help='full path to input file')
parser.add_argument('-m', '--mooring', nargs='+', 
                    help='add mooring location, individual moorings or [all] for all moorings \
                    must be after the infile if listing individual moorings [ie -m c1 c2 m2 m4]')
parser.add_argument('-ex', '--extents', nargs=1,
                    help='chooses extents of map, options are bering, chukchi, custom, and default')


args=parser.parse_args()

moorings = {'c1': [70.838,163.125], 'c2': [71.231,164.223], 'c3': [71.828,166.070],
            'c4': [71.038,160.514], 'c5': [71.203,158.011], 'c12': [67.911,168.195],
            'm2': [56.869,164.050], 'm4': [57.895,168.878], 'm5': [59.911,171.731],
            'm8': [62.194,174.688], 'm14': [64.0,167.929], 'k1': [57.855,163.667],
            'k2': [58.262,163.495], 'k3': [58.785,163.278], 'c9': [72.464,156.548]}

#extent note that easter extent is + deg E long
extent = {'chukchi': [180,210,64,73], 'bering': [174,206,51,66], 
          'default': [165,220,50,75],'custom':[185,210,56,61]} #adjust custom to desired extents

if args.extents:
    map_extents=extent[args.extents[0]]
else:
    map_extents=extent['default']

#get the date of the filename for future use
match = re.search(r'\d{8}', args.infile)
date = datetime.strptime(match.group(0), "%Y%m%d")
str_date = date.strftime("%Y-%m-%d")



ds = xr.open_dataset(args.infile)

ice_conc = ds['ICECON'].squeeze()

# 2. Define the projection (EPSG:3411)
proj_params = Proj("epsg:3411")

# 3. Create the 2D meshgrid of (x, y) coordinates
X_mesh, Y_mesh = np.meshgrid(ds.x.values, ds.y.values)

# 4. Perform the inverse transformation (x, y to lon, lat)
lon_2d, lat_2d = proj_params(X_mesh, Y_mesh, inverse=True)

# 5. Add the new 2D coordinates to the xarray Dataset
# Note: You should use the original dimension names ('y', 'x') for mapping
ds = ds.assign_coords(
    {'latitude': (('y', 'x'), lat_2d), 
     'longitude': (('y', 'x'), lon_2d)}
)

# Define the bounds
lon_w, lon_e, lat_s, lat_n = map_extents

#convert lon_e to negative longitude for the mask
lon_e = lon_e - 180 - 180

# Create the boolean mask. 
# xarray automatically handles the broadcast across the (y, x) grid.

if lon_w > lon_e:
    lon_mask = (lon_2d >= lon_w) | (lon_2d <= lon_e)
if lon_w < lon_e:
    lon_mask = (lon_2d >= lon_w) & (lon_2d <= lon_e)

mask = lon_mask & (lat_2d >= lat_s) & (lat_2d <= lat_n)

# Apply the mask to the dataset.
# The .where() method sets values outside the mask to NaN, but keeps the original grid structure.
#ds_sliced = ds.where(mask, drop=True)

# Apply the mask to the ICECON data and drop rows/columns where ALL values are NaN
# This reduces the grid size to just what is needed for the plot.
#ice_conc_plot = ds['ICECON'].where(mask, drop=True)

# Convert to 0-100% and replace fill values (e.g., 255) with NaN for plotting
ice_conc_plot = ice_conc * 100
#ice_conc_plot = np.where((mask) & (ice_conc_plot <= 100), ice_conc_plot, np.nan)
ice_conc_plot = np.where((ice_conc_plot <= 100), ice_conc_plot, np.nan)
#now set all zero values to nan so it makes the plot look better
#that is the zeros won't get plotted as dark in the colorbar
ice_conc_plot = np.where(ice_conc_plot == 0, np.nan, ice_conc_plot)
# The drop=True argument removes any rows (y) or columns (x) 
# where ALL values are masked (NaN), significantly reducing the data size.


#this was used for plotting the data as points, not used in grid plot
flat_lats = ds.latitude.values.flatten()
flat_lons = ds.longitude.values.flatten()
ice = ds.ICECON.values.flatten()


#this was used for plotting the data as points, not used in pcolormesh plot
# df_ice = pd.DataFrame({'latitude':flat_lats, 'longitude':flat_lons, "ice_conc":ice})
# df_ice['ice_conc'] = df_ice['ice_conc'] * 100
# df_ice = df_ice.loc[df_ice['ice_conc'] <= 100]


clon = map_extents[0] + (map_extents[1] - map_extents[0])/2
clat = map_extents[3] - (map_extents[3] - map_extents[2])/2
sparallel1 = 55
sparallel2 = 70

albers_proj = ccrs.AlbersEqualArea(
    central_longitude=clon,
    central_latitude=clat,
    standard_parallels=(sparallel1, sparallel2)
)

fig = plt.figure(figsize=(10,8))
ax = fig.add_subplot(1, 1, 1, projection=albers_proj)


ax.set_extent(map_extents, crs=ccrs.PlateCarree())

ax.add_feature(cfeature.OCEAN, color='azure')
ax.add_feature(cfeature.LAND, color='beige')
# 2. Use xarray's built-in plotting function for pcolormesh
# We use lon_2d and lat_2d as the coordinates for the mesh
# We use the sliced data for the color values
# Note: pcolormesh automatically handles the transformation
p = ax.pcolormesh(
    lon_2d,            # X coordinates (M, N shape)
    lat_2d,            # Y coordinates (M, N shape)
    ice_conc_plot,     # C data values (M, N shape, with NaN outside the map_extents)
    cmap=cmocean.cm.ice, 
    vmin=0, 
    vmax=100, 
    shading='auto',    # CRITICAL: Fixes the dimension mismatch error
    transform=ccrs.PlateCarree() 
)
plt.colorbar(p, pad=.1)
ax.add_feature(cfeature.COASTLINE)
gl = ax.gridlines()
gl.xlocator = plticker.MultipleLocator(10)
gl.ylocator = plticker.MultipleLocator(5)
gl.bottom_labels = True
gl.left_labels = True

transformation=ccrs.PlateCarree()
if args.mooring:    
    if args.mooring[0] == 'all':
        selected_moorings = list(moorings.keys())
    else:
        selected_moorings = args.mooring
    for i in selected_moorings:
        lat = moorings[i][0]
        lon = 360-moorings[i][1]
        ax.plot(lon, lat, 'r.', transform=transformation)
        #plot_title = key.swapcase() + ' ' + plot_title
        label = i.swapcase()
        label_lat = lat
        label_lon = lon
        ax.text(label_lon, label_lat, label, horizontalalignment='right', 
                verticalalignment='bottom', transform=transformation)

#ax.scatter(df_ice.longitude, df_ice.latitude, color='red', marker='o', s=1, transform=ccrs.PlateCarree())
plot_title = "AMSR2 25km Ice Concentration on " + str_date
ax.set_title(plot_title, fontsize=16, fontweight='bold', pad=20)
output_filename = 'amsr2_ice_25km_' + date.strftime("%Y_%m_%d") + ".png"
plt.savefig(output_filename, dpi=300)
print(f"Plot saved to '{output_filename}'")

