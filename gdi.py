# reading files
import pygrib                        # Provides a high-level interface to the ECWMF ECCODES C library for reading GRIB files
import xarray as xr                  # Provides N-dimensional labeled arrays and datasets in Python

# for plotting
# for map
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

#for visual plotting
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors

# metpy to calculate atmo stuff
import metpy
import metpy.calc as mpcalc
import metpy.constants as mpconsts
from metpy.plots import colortables

# misc os and other
import math
import sys
import os
import datetime

# function to get var based on pressure lev and var
# must know level of interest and var

def get_inputdat_GDI(var, extent):

  var_950 = grib.select(name= var, typeOfLevel = 'isobaricInhPa', level = 950)[0]
  var_850 = grib.select(name= var, typeOfLevel = 'isobaricInhPa', level = 850)[0]
  var_700 = grib.select(name= var, typeOfLevel = 'isobaricInhPa', level = 700)[0]
  var_500 = grib.select(name= var, typeOfLevel = 'isobaricInhPa', level = 500)[0]

  var_950_ext = domain_ext(var_950, extent)[0]
  var_850_ext = domain_ext(var_850, extent)[0]
  var_700_ext = domain_ext(var_700, extent)[0]
  var_500_ext = domain_ext(var_500, extent)[0]
  return(var_950_ext, var_850_ext, var_700_ext, var_500_ext)

# sets a domain for variables

def domain_ext(var, extent):

  # Read the data for a specific region
  var_ext = var.data(lat1=extent[1],lat2=extent[3],lon1=extent[0]+360,lon2=extent[2]+360)

  return(var_ext)

# sets a domain for variables

def domain_ext(var, extent):

  # Read the data for a specific region
  var_ext = var.data(lat1=extent[1],lat2=extent[3],lon1=extent[0]+360,lon2=extent[2]+360)

  return(var_ext)

from metpy.units import units

# Image path
path = ('C:\Users\michell.tinoco\Desktop\GFS_FILES_EX/gfs.t00z.pgrb2.1p00.f006')

# open the GRIB file
grib = pygrib.open(path)

# Select the extent [min. lon, min. lat, max. lon, max. lat]
extent = [-110.0, 5, -55.00, 30.00]

# ============ INPUT DATA FOR ALGORITHM ============ #
# select var names
var_1 = 'Temperature'
var_2 = 'Relative humidity'

#use function get_dat to retrieve data @ appropriate hPa levs

# Read the temperature in 950 hPa
t_950, t_850, t_700, t_500 = get_inputdat_GDI(var_1, extent)*units.kelvin
print('temperature data read succesfully!')

#convert to celcius (does work)
#t_950 = t_950.to('degC')

# Read the relative humidity in 950 hPa
rh_950, rh_850, rh_700, rh_500 = get_inputdat_GDI(var_2, extent)

# surface pressure
sfc_pres = grib.select(name='Surface pressure')[0]
sfc_pres_ext, lats, lons = sfc_pres.data(lat1=extent[1],lat2=extent[3],lon1=extent[0]+360,lon2=extent[2]+360)

#convert pa to hpa
sfc_pres_ext = sfc_pres_ext*units.Pa
psfc = sfc_pres_ext.to('hPa')


print('relative humidity data read succesfully!')

from metpy.units import units
from metpy.calc import mixing_ratio_from_relative_humidity
import numpy as np

# ============ Layers ============ #
# Constanst for layers A, B, & C
alpha = units.Quantity(-10 , 'K') #limits excessive GDI in regions w/ lots of moisture above 850hPa
L0 = units.Quantity(2.69 * 10**6 , ' J kg^-1') # latent heat of vapor.
Cpd = units.Quantity(1005.7 , 'J kg^-1 K-1') # specific heat of dry air

#equivalent potential temperature temp
pres = [950, 850, 700, 500]
temps = [t_950, t_850, t_700, t_500]
relhum = [rh_950, rh_850, rh_700, rh_500]

theta = []
mixrat = []
for i in range(len(pres)):
  #calculating theta
  thet = metpy.calc.potential_temperature(pres[i]* units.hPa, temps[i])
  theta.append(thet)

  #calculating mixing ration
  mxrt = mixing_ratio_from_relative_humidity(pres[i]* units.hPa, temps[i], relhum[i]* units.percent).to('g/kg')
  mixrat.append(mxrt)

ept_a = theta[0] * np.exp((L0 * mixrat[0])/(Cpd * t_850))
ept_b = ((0.5 * (theta[1] + theta[2])) * np.exp((L0 * (0.5 * (mixrat[1] + mixrat[2]))) / (Cpd * t_850))) + alpha
ept_c = (theta[3] * np.exp((L0 * mixrat[3]) / (Cpd * t_850)) ) + alpha

print(' Equivalent Potential Temp Calculations Successfully Calculated! ')

# CBI Calculation
#constants
gamma = units.Quantity(6.5*10**-2 , 'K^-2') #scaling constant. #tentative change but k^-1
beta = units.Quantity(303 , 'K') #why 303k?

#mid level
ME = ept_c - beta

#low level
LE = ept_a - beta

#conditional
CBI = np.where(LE <= 0, 0, gamma * LE * ME)
print('CBI Units:', (LE * ME).units)  # Print the units of CBI
print(' CBI successfully calculated! ')


# MWI
#constants
mu = units.Quantity(-7 , 'K^-1') #scaling constant to set MWI as neg
tau = units.Quantity(263.15 , 'K') #temp denoting 500hPa temps (but why only -10C)

#conditional
MWI = np.where((t_500 - tau) <= 0, 0, mu * (t_500 - tau))
print('MWI Units:', MWI.units)  # Print the units of MWI
print(' MWI successfully calculated! ')

# II

#constants
sigma = units.Quantity(1.5 , 'K^-1') #scaling constant

#calculations
s = t_950 - t_700
d = ept_b - ept_a

#conditional
II = np.where((s + d) <= 0, sigma * (s + d), 0)
print('II Units:', II.units)  # Print the units of II
print(' II successfully calculated! ')


#TC

# constants
p1 = units.Quantity(500 , 'hPa')
p2 = units.Quantity(9000 , 'hPa')
p3 = 18

#calculations
TC = p3 - (p2 / ( psfc - p1))
print('TC Units:', TC.units)  # Print the units of TC
print(' TC successfully calculated! ')

# GDI

GDI = CBI + MWI + II + TC
print('GDI Units:', GDI.units)  # Print the units of GDI
print(' GDI successfully calculated! ')


#Information for plotting

#uses sfc_pres to get this information. can be any read data tbh
# Get information from the file
init  = str(sfc_pres.analDate)      # Init date / time
run   = str(sfc_pres.hour).zfill(2) # Run
ftime = str(sfc_pres.forecastTime)  # Forecast hour
valid = str(sfc_pres.validDate)     # Valid date / time
print('Init: ' + init + ' UTC')
print('Run: ' + run + 'Z')
print('Forecast: +' + ftime)
print('Valid: ' + valid + ' UTC')

#plotting

# Create a figure and axis using Cartopy's PlateCarree projection
fig, ax = plt.subplots(figsize = (10, 10), subplot_kw={'projection': ccrs.PlateCarree()})

# Select the extent [min. lon, min. lat, max. lon, max. lat]
extent_map = [extent[0], extent[2], extent[1], extent[3]]

# Add features to the map
ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth = 0.8, color = 'black')  # Make the borders semi-transparent
ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth = 0.8, color = 'black')  # Make the borders semi-transparent

# Add gridlines
gl = ax.gridlines(crs=ccrs.PlateCarree(), color='gray', alpha=1.0, linestyle='--', linewidth=0.25, xlocs=np.arange(-180, 180, 5), ylocs=np.arange(-90, 90, 5), draw_labels=True)
gl.top_labels = False
gl.right_labels = False

# Add a title
plt.title(f'GFS: Galvez Davison Index (GDI): {init}' , fontweight='bold', fontsize=10, loc='left')
plt.title('Valid: ' + valid  , fontsize=10, loc='right')

#colorbar
# Define de contour interval
data_min = -30
data_max = 70
interval = 5
levels = np.arange(data_min,data_max,interval)

# Create the color scale
colors = ["#323232", "#646464", "#737373", "#7e7e7e", "#909090", "#a3a3a3", "#b1b1b1", "#bcbcbc", "#bbc7cb", "#b2d2dd", "#90d5bb", "#55d065", "#5acf28", "#bad411", "#ffcc00", "#ffa900", "#fc8106", "#eb4722", "#d8133a", "#ac0a1d"]
cmap = matplotlib.colors.ListedColormap(colors)
cmap.set_over('#800000')
cmap.set_under('#000000')

# ============ data ============ #

img1 = ax.contourf(lons, lats, GDI, extent = extent_map, cmap = cmap, levels = levels)
img2 = ax.contour(lons, lats, GDI, colors='white', linewidths=0.3, levels=[-5, 5, 15, 25, 35, 45, 55, 65])
#ax.clabel(img2, inline=1, inline_spacing=0, fontsize='10',fmt = '%1.0f', colors= 'white')

#colorbar
from mpl_toolkits.axes_grid1 import make_axes_locatable

divider = make_axes_locatable(ax)
# Pass the projection to the append_axes function
cax = divider.append_axes('bottom', size='5%', pad=0.05, axes_class=plt.Axes) # Pass axes_class=plt.Axes
cbar = fig.colorbar(img1, cax=cax, orientation='horizontal', extend='both')

plt.savefig('GDI_example.png')
