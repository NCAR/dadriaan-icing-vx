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
CLIPMIN = -0.1 
CLIPMAX = 1.1
#DSTRING = "20190201"
#FSTRING = "20190201_0008*.nc"
DSTRING = sys.argv[1]
FSTRING = sys.argv[2]
FIELD = 'ICE_PROB'
UNITS = '%'
BINS = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2]
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

# If we want to convert probability to potential, do it here
if PROB2POT:
  ncA[FIELD] = ncA[FIELD]/0.85
  ncB[FIELD] = ncB[FIELD]/0.85

# Replace missing values if requested
if subMissing:
  varA = ncA[FIELD].fillna(subMissingA)
  varB = ncB[FIELD].fillna(subMissingB)

  # Make sure there are no fill values, there shouldn't be
  #missValsA = ncA.where(xr.ufuncs.isnan(difference),drop=True)
  #missValsB = ncB.where(xr.ufuncs.isnan(difference),drop=True)

# Clip values prior to histogram
varAClip = varA.clip(min=CLIPMIN,max=CLIPMAX)
varBClip = varB.clip(min=CLIPMIN,max=CLIPMAX)

# Total number of any points
nPtsVol = int(np.prod(list(ncA.sizes.values())))

# Set a new figure window
#fig = plt.figure(1, figsize=(22, 15))
fig = plt.figure(1, figsize=(15, 15))

# Plot
ax1 = plt.subplot(111)
#ax1.set_ylim((0,1.0e7))

# Histogram
Acounts, Abins, Abars = xr.plot.hist(varAClip,ax=None,bins=BINS)
Bcounts, Bbins, Bbars = xr.plot.hist(varBClip,ax=None,bins=BINS)

# Remove histogram bars from figure
At = [b.remove() for b in Abars]
Bt = [b.remove() for b in Bbars]

# Plot values with bar
barA = Acounts
barB = Bcounts
xA = np.arange(len(barA))
xB = xA+0.25
if DEBUG:
  print("XA TICKS:")
  print(xA)
  print("XB TICKS:")
  print(xB)

# Set xtick locations and labels
xtick = []
barWidth = 0.5
groupWidth = 1.5 # 
tickOff = 0.1 # Offset from edge of bar group
count = 0
bins = BINS
for v in xA:

  # The current v
  if DEBUG:
    print(v)
  
  # The current bin
  if DEBUG:
    print(bins[v])

  # The tick left of the current bar group
  if DEBUG:
    print(v-(barWidth/2.0)-tickOff)
  xtick.append(v-(barWidth/2.0)-tickOff)

  # The tick right of the current bar group
  #print((v-(barWidth)/2.0)+groupWidth+tickOff)
  #xtick.append((v+0.25-(barWidth)/2.0)+groupWidth+tickOff)
  #print(v-(barWidth)/2.0)-tickOff+groupWidth+tickOff)
  #xtick.append

  v = v + 1

# Append one last tick mark
xtick.append(max(xtick)+1.0)
  
ax1.set_xticks(xtick)
ax1.set_xticklabels(bins)
ax1.set_xlabel('%s (%s)' % (FIELD,UNITS))
ax1.set_ylabel('Count')

if DEBUG:
  print(barA)
  print(barB)
  print(np.abs(barA-barB))

# Plot A over B, so plot B first, then A
bpB = plt.bar(xB,barB,width=0.50,color='black',edgecolor='black',linewidth=0.5)
bpA = plt.bar(xA,barA,width=0.50,color='cyan',edgecolor='black',linewidth=0.5)
#bpA = plt.bar(xA,barA-barB,width=0.50)
ax1.relim()

# Compute bin percentages, compared to total points with any differences
#binPercs = (hist[0]/nPtsDiff)*100.0

# Add labels to the top of each bar
#binCnt = 0
#for p in hist[2]:
#  txX = p.get_x()+p.get_width()/2.0
#  txY = p.get_height()+0.05*max(binPercs)
#  ax1.text(txX,txY,'%0.2f' % (binPercs[binCnt]),ha='center')
#  binCnt = binCnt + 1

# Figure title
#fig.suptitle('TESTING')
ax1.set_title('Histogram of dataset A (cyan) and dataset B (black)\n Nfiles: %d Nsample: %d' % (len(Afiles),int(nPtsVol)))

# Save the image
plt.savefig('compare.png')
