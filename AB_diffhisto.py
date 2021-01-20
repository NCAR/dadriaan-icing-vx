#!/usr/bin/env python

# Imports
import xarray as xr
import glob, time, datetime, sys
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import numpy as np

# Vars
#chunks = {'y0':158,'x0':158}
chunks = {'time':1,'z0':1}
subMissing = True
subMissingA = 10000.0
subMissingB = 10000.0
DIFFTHRESH = 0.1 # Threshold for considering the values different between versions
DIFFOK = 5 # Number of differences to allow at any i,j over the analysis period (i.e. mask out anything <= this value)
DEBUG = False
DIFFDIR = "amb" # or "bma"
DIFFTYPE = "all" # or "pos" or "neg"
#DSTRING = "20190201" # YYYYMMDD or "*" wildcard for all
DSTRING = sys.argv[1]
#FSTRING = "*.nc"
FSTRING = sys.argv[2]
CLIPMIN = -1.1
CLIPMAX = 1.1
FIELD = 'ICE_PROB'
UNITS = '%'
BINS = [-1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2]
TICKS = [-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1]
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

# Compute the total number of points being analyzed
nptsanly = float(np.prod(list(ncB.sizes.values())))

# If we want to convert probability to potential, do it here
if PROB2POT:
  ncA[FIELD] = ncA[FIELD]/0.85
  ncB[FIELD] = ncB[FIELD]/0.85

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

# Clip the difference values
diffClip = difference.clip(min=CLIPMIN,max=CLIPMAX)
#print(diffClip.max().values)

# Mask where differences are >= some thresh. Drop any points where they don't meet the thresh.
print("\nMASKING DIFFERENCE...")
if DIFFTYPE=='all':
  diffmask = diffClip.where(((diffClip>=DIFFTHRESH) | (diffClip<=-1.0*DIFFTHRESH)),drop=True)
elif DIFFTYPE=='pos':
  diffmask = diffClip.where(diffClip>=DIFFTHRESH,drop=True)
elif DIFFTYPE=='neg':
  diffmask = diffClip.where(diffClip<=-1.0*DIFFTHRESH,drop=True)
else:
  print("\nFATAL! UNKNOWN DIFFTYPE")
  exit()

# Total number of points with any difference
# First compute the total number of points remaining in the masked array
nDiff = int(np.prod(list(diffmask.sizes.values())))

# Now subtract off the total number of points that are NaN. These are points we don't want
nPtsDiff = nDiff - xr.ufuncs.isnan(diffmask).sum().values

# Total number of any points
nPtsVol = int(np.prod(list(diffClip.sizes.values())))

# Set a new figure window
#fig = plt.figure(1, figsize=(22, 15))
fig = plt.figure(1, figsize=(15, 15))

# Plot
ax1 = plt.subplot(111)
ax1.set_xticks(ticks=TICKS)

# Histogram
hcounts, hbins, hbars = xr.plot.hist(diffmask,ax=ax1,bins=BINS,rwidth=0.75)

ax1.set_xlabel('%s (%s)' % (FIELD,UNITS))
ax1.set_ylabel('Count')

# Compute bin percentages, compared to total points with any differences
binPercs = (hcounts/nPtsDiff)*100.0

# Add labels to the top of each bar
binCnt = 0
for p in hbars:
  txX = p.get_x()+p.get_width()/2.0
  txY = p.get_height()+(0.01*max(hcounts))
  ax1.text(txX,txY,'%0.2f' % (binPercs[binCnt]),ha='center')
  binCnt = binCnt + 1

# Figure title
#fig.suptitle('Differences > %03f %%, N = %15d' % (DIFFTHRESH,int(np.prod(list(diffmask.sizes.values())))))
ax1.set_title('Differences > %f %s, nHrs = %04d\n%f %% of all points' % (DIFFTHRESH,UNITS,int(len(Afiles)),(nPtsDiff/nptsanly)*100.0))

# Save the image
plt.savefig('diffhisto.png')
