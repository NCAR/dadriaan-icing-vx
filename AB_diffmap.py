#!/usr/bin/env python

# Imports
import xarray as xr
import glob, time, datetime, sys
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import BoundaryNorm
from matplotlib.colors import Normalize
#import matplotlib as mpl
#mpl.rcParams['figure.dpi'] = 1000

# Vars
#chunks = {'y0':158,'x0':158}
chunks = {'time':1,'z0':1}
subMissing = True
subMissingA = 10000.0
subMissingB = 10000.0
DIFFTHRESH=0.1 # Threshold for considering the values different between versions
DIFFOK = 5 # Number of differences to allow at any i,j over the analysis period (i.e. mask out anything <= this value)
DEBUG = False
DIFFDIR = "amb" # or "bma"
DIFFTYPE = "all" # or "pos" or "neg"
#DSTRING = "20190202" # YYYYMMDD or "*" wildcard for all
DSTRING = sys.argv[1]
#FSTRING = "*.nc"
FSTRING = sys.argv[2]
#LEVELS = np.arange(0,320,10)
LEVELS = np.array([0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0])*100.0
#LEVELS = None
MASKTYPE = "miss" # numeric value (i.e. 0) or "miss"
BORDERCOLOR = 'black'
BORDERWIDTH = 0.1
#MYDPI = 300
#MYDPI = 96
FIELD = 'ICE_PROB'
UNITS = '%s'
PROB2POT = True

# File paths
# FIRST
#Apath = "/d1/dadriaan/data/versionB/netcdf/conus_ICE_PROB" # BASELINE
#Bpath = "/d1/dadriaan/data/versionC/ral-ifi-cip/main/data_dir/num25_hrrr/data/netcdf/cip/pressure/conus_ICE_PROB" # NOBLEND

# SECOND
#Apath = "/d1/dadriaan/data/versionC/ral-ifi-cip/main/data_dir/num25_hrrr/data/netcdf/cip/pressure/conus_ICE_PROB" # NOBLEND
#Bpath = "/var/autofs/mnt/khaba_d1/dadriaan/data/versionC/ral-ifi-cip/main/data_dir/num25_hrrr/data/netcdf/cip/pressure/conus_ICE_PROB" # RADIUS

# THIRD
Apath = "/var/autofs/mnt/khaba_d1/dadriaan/data/versionC/ral-ifi-cip/main/data_dir/num25_hrrr/data/netcdf/cip/pressure/conus_ICE_PROB" # RADIUS
Bpath = "/var/autofs/mnt/khaba_d1/dadriaan/data/versionC/ral-ifi-research/glm_test/data_dir/num25_hrrr/data/netcdf/cip/pressure/conus_ICE_PROB" # GLM

# Collect the files
Afiles = glob.glob('%s/%s/%s' % (Apath,DSTRING,FSTRING))
Bfiles = glob.glob('%s/%s/%s' % (Bpath,DSTRING,FSTRING))
#print(Afiles)
#print(Bfiles)
#Afiles.sort()
#Bfiles.sort()
#for a in Afiles:
#  print(a)
#for b in Bfiles:
#  print(b)

# Make sure there's the same number of files for each version
if len(Afiles) != len(Bfiles):
  print("\nFATAL! Number of files different.")
  print("nfA = %04d" % (len(Afiles)))
  print("nfB = %04d" % (len(Bfiles)))
  print("\n")
  exit()

# Open the files with Xarray
t1 = datetime.datetime.now()
print("READING A...")
#ncA = xr.open_mfdataset(Afiles,combine='by_coords')
ncA = xr.open_mfdataset(Afiles,chunks=chunks,combine='by_coords')
t2 = datetime.datetime.now()
print(t2-t1)
print(ncA)
print("\nREADING B...")
#ncB = xr.open_mfdataset(Bfiles,combine='by_coords')
ncB = xr.open_mfdataset(Bfiles,chunks=chunks,combine='by_coords')
t3 = datetime.datetime.now()
print(t3-t2)
print(ncB)
#print(ncA['grid_mapping_0'])
#print(ncB['grid_mapping_0'])

# If we want to convert probability to potential, do it here
if PROB2POT:
  ncA[FIELD] = ncA[FIELD]/0.85
  ncB[FIELD] = ncB[FIELD]/0.85

# Figure out maxdiffs
maxdiff = ncA.dims['time']*ncA.dims['z0']

# Take a difference
print("\nCOMPUTING DIFFERENCE...")
if DIFFDIR=='amb':
  if subMissing:
    difference = ncA[FIELD].fillna(subMissingA)-ncB[FIELD].fillna(subMissingB)
  else:
    difference = ncA[FIELD]-ncB[FIELD]
else:
  if subMissing:
    difference = ncB[FIELD].fillna(subMissingB)-ncA[FIELD].fillna(subMissingA)
  else:
    difference = ncB[FIELD]-ncA[FIELD]

# Make sure there are no fill values, there shouldn't be
#missVals = difference.where(xr.ufuncs.isnan(difference),drop=True)

# First, we want to mask where differences are >= some thresh
print("\nMASKING DIFFERENCE...")
if DIFFTYPE=='all':
  diffmask = xr.where(((difference>=DIFFTHRESH) | (difference<=-1.0*DIFFTHRESH)),1.0,0.0)
elif DIFFTYPE=='pos':
  diffmask = xr.where(difference>=DIFFTHRESH,1.0,0.0)
elif DIFFTYPE=='neg':
  diffmask = xr.where(difference<=-1.0*DIFFTHRESH,1.0,0.0)
else:
  print("\nFATAL! UNKNOWN DIFFTYPE")
  exit()

# Condense differences into a 2D count
print("\nSUMMARIZING DIFFERENCES...")
diffsum = diffmask.sum(dim='time').sum(dim='z0')

print("\nMASKING ZEROES...")
if MASKTYPE=='zero':
  # Set to value
  diffsum = xr.where(diffsum>DIFFOK,diffsum,0)
elif MASKTYPE=='miss':
  # Set to 'NA'
  diffsum = diffsum.where(diffsum>DIFFOK)
else:
  print("\nFATAL! UNKNOWN MASKTYPE")
  exit()

# Four tests
# 1. diffsum with NA Lambert and Latlon
# 2. diffsum with 0 Lambert and Latlon
# 3. ICE_PROB with Lambert and Latlon
# 4. same as 2, but control colors shown with contour levels?

# Test group 1
# mask with NA, Lambert pcolor
# mask with NA, Latlon pcolor
# mask with NA, Lambert contourf --> this fails
# mask with NA, Latlon contourf

# Test group 2
# mask with 0, Lambert pcolor
# mask with 0, Latlon pcolor
# mask with 0, Lambert contourf
# mask with 0, Latlon contourf

# Test group 3
# prob, Lambert pcolor
# prob, Latlon pcolor
# prob, Lambert contourf
# prob, Latlon contourf

# Test group 4
# mask with NA, Latlon contourf, levels
# mask with 0, Lambert contourf, levels
# mask with 0, Latlon contourf, levels

# Testing results...we probably want to mask with zeros, and then exclude them at the contouring stage via
#                   the colormap or by masking them out after the fact. I think including missing values
#                   affects the contouring algorithm, and since the values aren't actually missing we should
#                   be including them for the contouring algorithm to use. If data truly ARE missing, we will
#                   have to test what to do like for icing output fields.
#
# Testing results...we also probably just want to use a lat/lon projection since it seems defining the projection
#                   as LambertConformal then translating back to PlateCarree() takes forever.
#
# Testing results...pcolormesh produces similar results to contourf but slightly muted. This may actually be correct,
#                   since contourf bins the results and pcolormesh interpolates to a very smooth grid making it
#                   visually different from a color perspective but yet potentially accurate from a numerical standpoint.
#                   pcolormesh is loads faster, so perhaps we should investigate whether this is OK to use.

# Coordinate reference system
print("\nCREATING PLOT...")
#crs = ccrs.LambertConformal(central_longitude=-95.0, central_latitude=25.0) # NCL coords, old RAP grid
crs = ccrs.LambertConformal(central_longitude=-97.5, central_latitude=38.5) # HRRR coords
#crs = ccrs.PlateCarree() # Lat/Lon

# Turn off showing plot?
#plt.ioff()

# Set a new figure window in inches
#fig = plt.figure(1, figsize=(22, 15))
#fig = plt.figure(1, figsize=(15, 15))
fig = plt.figure(1, figsize=(19.5, 19.5))
#fig = plt.figure(1, figsize=(6.4, 4.8))
#fig = plt.figure(1, figsize=(100, 100))
#fig = plt.figure(1)

# Set a new figure window in pixels
#px = 1/plt.rcParams['figure.dpi']
#px = 1/MYDPI
#print(plt.rcParams['figure.dpi'])
#print(px)
#print(1200*px)
#fig = plt.figure(1, figsize=(1200*px,1200*px))
#fig = plt.figure(1, figsize=(2000*px,2000*px))

# Figure attributes
#print(fig.get_size_inches())
#print(fig.dpi)
#print(fig.get_window_extent())

# Setup axis
def axis_setup(ax):
  # Zoom the map
  #ax.set_extent([235.,290.,20.,55.],ccrs.PlateCarree())
  #ax.set_extent([230.,300.,16.,60.],ccrs.PlateCarree())
  #ax.set_extent([235.,290.,20.,55.])
  #ax.set_global()
  #ax.set_extent([270.0,280.0,40.0,45.0]) # TEST ZOOM Lake Erie/MI
  #ax.set_extent([250.0,260.0,30.0,35.0]) # TEST ZOOM AZ

  # Add borders and map features
  ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=BORDERWIDTH,edgecolor=BORDERCOLOR)
  ax.add_feature(cfeature.STATES, linewidth=BORDERWIDTH,edgecolor=BORDERCOLOR)
  ax.add_feature(cfeature.BORDERS, linewidth=BORDERWIDTH,edgecolor=BORDERCOLOR)
  #ax.coastlines()
  return ax

# Get some more map outlines to use

# Get a slice of data to plot
#plotData = ncA.ICE_PROB.isel(time=0,z0=0)

# Plot setup
#ax1 = plt.subplot(111,projection=crs)
ax1 = plt.axes([0,0,1,1],projection=crs)
#ax1 = plt.axes([0,0,1,1],projection=crs)
axis_setup(ax1)

# Stop for debugging
#import pdb; pdb.set_trace()

# Difference plotting
# Raw counts divided by total max differences possible 
#norm = Normalize(0.0,1.0,clip=True) # Normalize 0 to 1.0 if we're dividing by maxdiff
norm = BoundaryNorm(LEVELS,ncolors=plt.get_cmap('plasma').N,clip=True)
cf1 = ax1.pcolormesh(diffsum.lon0,diffsum.lat0,(diffsum/maxdiff)*100.0,cmap='plasma',transform=ccrs.PlateCarree(),shading='nearest',norm=norm)

# Raw counts
#norm = Normalize(0.0,maxdiff,clip=True) # Normalize 0 to maxdiff?
#cf1 = ax1.pcolormesh(diffsum.lon0,diffsum.lat0,diffsum,cmap='plasma',transform=ccrs.PlateCarree(),shading='nearest')

# CONTOURF TEST
#cf1 = ax1.contourf(diffsum.lon0,diffsum.lat0,diffsum,cmap='plasma',transform=ccrs.PlateCarree(),levels=LEVELS)

# Data plotting
#cf1 = ax1.contourf(plotData.lon0,plotData.lat0,plotData,cmap='plasma',transform=ccrs.PlateCarree())
#cf1 = ax1.pcolormesh(plotData.lon0,plotData.lat0,plotData,cmap='plasma',transform=ccrs.PlateCarree())

ax1.set_title('Max possible differences at each grid cell: %d' % (int(maxdiff)),loc='right')
ax1.set_title('Geographic differences of %s in height and time for %d files' % (FIELD,len(Afiles)),loc='left')

cb1 = fig.colorbar(cf1,ax=ax1,orientation='horizontal',shrink=0.74,pad=0.01,ticks=LEVELS)
#cb1 = fig.colorbar(cf1,ax=ax1,orientation='horizontal',shrink=0.74,pad=0.01)
#cb1.set_xticks([0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])
#cb1.set_xticklabels(
cb1.set_label('%% of total points with %s differences > %f' % (FIELD,DIFFTHRESH),size='x-large')
#fig.suptitle('DIFFERENCES',y=0.90)
#plt.show()
#plt.savefig('testing%s.png' % (DSTRING),dpi=MYDPI)
plt.savefig('testing%s.png' % (DSTRING),bbox_inches='tight')
