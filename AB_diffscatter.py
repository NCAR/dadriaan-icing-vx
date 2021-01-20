#!/usr/bin/env python

# Imports
import xarray as xr
import glob, time, datetime, sys
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import numpy as np
#import seaborn as sns

np.set_printoptions(suppress=True)

# Set boolens to request plots
# B2 = Total occurrences of A and B both equal to 0.0
# B3 = Total occurrences of A changing from missing to 0.0 in B
# B4 = Histogram of values that were > 0.0 in A but 0.0 in B
# C2 = Total occurrences of A changing from 0.0 to missing in B
# C3 = Total occurrences of A and B both equal to missing
# C4 = Histogram of values that were > 0.0 in A that became missing in B
# C5 = Histogram of values that were valid in A (>=0.0) but missing in B
# D2 = Histogram of values that were 0.0 in A but > 0.0 in B
# D3 = Histogram of values that were missing in A but > 0.0 in B
# D4 = X-Y plot of any values > 0.0 in both A and B
# E3 = Histogram of values that were missing in A, but are now valid in B (>=0.0)
# E5 = X-Y plot of any values >= 0.0 in both A and B
B2=B3=B4=C2=C3=C4=C5=D2=D3=D4=E3=E5=True # ALL

# Vars
#chunks = {'y0':158,'x0':158}
chunks = {'time':1,'z0':1}
subMissing = False
subMissingA = 2.0
subMissingB = -2.0
DIFFTHRESH = 0.1 # Threshold for considering the values different between versions
DIFFOK = 5 # Number of differences to allow at any i,j over the analysis period (i.e. mask out anything <= this value)
DEBUG = False
DIFFDIR = "amb" # or "bma"
DIFFTYPE = "all" # or "pos" or "neg"
#CLIPMIN = -1.1
#CLIPMAX = 1.1
#CLIPMIN = -2.0
#CLIPMAX = 2.0
CLIPMIN = -0.85
CLIPMAX = 0.85
#DSTRING = "20190201"
#FSTRING = "20190201_0008*.nc"
DSTRING = sys.argv[1]
FSTRING = sys.argv[2]
#BINS = [0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9]
#BINS = np.arange(0.0,0.9,0.01)
BINS = np.arange(0.0,1.0,0.01)
DENS = True
YMIN = 0.0
YMAX = 100.0
FIGSIZE = (15,15)
FIELD = 'ICE_PROB'
UNITS = '%s'
PROB2POT = True

# File paths
# FIRST
#Apath = "/d1/dadriaan/data/versionB/netcdf/conus_ICE_PROB" # BASELINE
#Bpath = "/d1/dadriaan/data/versionC/ral-ifi-cip/main/data_dir/num25_hrrr/data/netcdf/cip/pressure/conus_ICE_PROB" # NOBLEND

# SECOND
Apath = "/d1/dadriaan/data/versionC/ral-ifi-cip/main/data_dir/num25_hrrr/data/netcdf/cip/pressure/conus_ICE_PROB" # NOBLEND
Bpath = "/var/autofs/mnt/khaba_d1/dadriaan/data/versionC/ral-ifi-cip/main/data_dir/num25_hrrr/data/netcdf/cip/pressure/conus_ICE_PROB" # RADIUS

# THIRD
#Apath = "/var/autofs/mnt/khaba_d1/dadriaan/data/versionC/ral-ifi-cip/main/data_dir/num25_hrrr/data/netcdf/cip/pressure/conus_ICE_PROB" # RADIUS
#Bpath = "/var/autofs/mnt/khaba_d1/dadriaan/data/versionC/ral-ifi-research/glm_test/data_dir/num25_hrrr/data/netcdf/cip/pressure/conus_ICE_PROB" # GLM

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

# **D4** Subset to only values > 0.0 in A and B
# Q: What does the scatter plot of positive A vs positive B look like (> 0.0 & > 0.0)
if D4:
  print("======D4======")
  fig1,ax1 = plt.subplots(figsize=FIGSIZE)
  tx = datetime.datetime.now()
  d4A = ncA[FIELD].where((ncA[FIELD]>0.0) & (ncB[FIELD]>0.0))
  print(d4A.count().values)
  print(datetime.datetime.now()-tx)
  d4B = ncB[FIELD].where((ncA[FIELD]>0.0) & (ncB[FIELD]>0.0))
  print(datetime.datetime.now()-tx)
  difftest = d4A-d4B
  diffA = d4A.where(((difftest>=DIFFTHRESH) | (difftest<=-1.0*DIFFTHRESH)),drop=True)
  npda = int(np.prod(list(diffA.sizes.values())))
  ndiffA = npda - xr.ufuncs.isnan(diffA).sum().values
  print(datetime.datetime.now()-tx)
  diffB = d4B.where(((difftest>=DIFFTHRESH) | (difftest<=-1.0*DIFFTHRESH)),drop=True)
  print(datetime.datetime.now()-tx)
  print(diffA)
  print(diffB)
  ax1.set_ylim(bottom=0.0,top=1.0)
  ax1.set_xlim(left=0.0,right=1.0)
  ax1.set_xlabel('B values')
  ax1.set_ylabel('A values')
  #ax1.scatter(diffA,diffB,marker='.',c=['red'])
  cnt, xed, yed, h2d1 = ax1.hist2d(diffB.values.flatten(),diffA.values.flatten(),bins=BINS,density=False)
  cb1 = fig1.colorbar(h2d1, ax=ax1, fraction=0.045, pad=0.05)
  cb1.set_label('Count')
  ax1.set_title('Values > 0.0 in A and B\nN = %d\n%f %% of all points' % (ndiffA,(ndiffA/nptsanly)*100.0))
  fig1.savefig('d4_xy.png')
  print(datetime.datetime.now()-tx)

# **E5** Count of valid in A and valid in B
# Q: What does the scatter plot of valid A vs valid B look like (>= 0.0 & >= 0.0)
if E5:
  print("======E5======")
  fig2,ax2 = plt.subplots(figsize=FIGSIZE)
  tx = datetime.datetime.now()
  e5A = ncA[FIELD].where((ncA[FIELD]>=0.0) & (ncB[FIELD]>=0.0))
  print(e5A.count().values)
  print(datetime.datetime.now()-tx)
  e5B = ncB[FIELD].where((ncB[FIELD]>=0.0) & (ncA[FIELD]>=0.0))
  print(datetime.datetime.now()-tx)
  difftest = e5A-e5B
  diffA = e5A.where(((difftest>=DIFFTHRESH) | (difftest<=-1.0*DIFFTHRESH)),drop=True)
  npda = int(np.prod(list(diffA.sizes.values())))
  ndiffA = npda - xr.ufuncs.isnan(diffA).sum().values
  print(datetime.datetime.now()-tx)
  diffB = e5B.where(((difftest>=DIFFTHRESH) | (difftest<=-1.0*DIFFTHRESH)),drop=True)
  print(datetime.datetime.now()-tx)
  ax2.set_ylim(bottom=0.0,top=1.0)
  ax2.set_xlim(left=0.0,right=1.0)
  ax2.set_xlabel('B values')
  ax2.set_ylabel('A values')
  #ax2.scatter(diffA.values,diffB.values,marker='.',c=['red'])
  cnt, xed, yed, h2d2 = ax2.hist2d(diffB.values.flatten(),diffA.values.flatten(),bins=BINS,density=False)
  cb2 = fig2.colorbar(h2d2, ax=ax2, fraction=0.045, pad=0.05)
  cb2.set_label('Count')
  ax2.set_title('Values >= 0.0 in A and B\nN = %d\n%f %% of all points' % (ndiffA,(ndiffA/nptsanly)*100.0))
  fig2.savefig('e5_xy.png')
  print(datetime.datetime.now()-tx)

# **C5** Histogram of valid A values where B is missing
# Q: If there are missing values in B that used to be valid in A, what is the distribution of those A values?
if C5:
  print("======C5======")
  fig3,ax3 = plt.subplots(figsize=FIGSIZE)
  tx = datetime.datetime.now()
  c5 = ncA[FIELD].where(((ncA[FIELD]>=0.0) & (xr.ufuncs.isnan(ncB[FIELD])))).values.flatten()
  nc5 = ((ncA[FIELD]>=0.0) & (xr.ufuncs.isnan(ncB[FIELD]))).sum().values
  print(datetime.datetime.now()-tx)
  ax3.set_ylim(bottom=YMIN,top=YMAX)
  ax3.set_title('A >= 0.0 where B is missing\nN = %d\n%f %% of all points' % (nc5,(nc5/nptsanly)*100.0))
  ax3.set_xlabel('Values')
  ax3.set_ylabel('Percent of all values meeting condition')
  c5h = ax3.hist(c5,bins=BINS,density=DENS)
  print(c5h[0])
  print(datetime.datetime.now()-tx)
  fig3.savefig('c5_hg.png')

# **E3** Histogram of valid B values where A is missing
# Q: If there were missing values in A that are now valid in B, what is the distribution of those B values?
if E3:
  print("======E3======")
  fig4,ax4 = plt.subplots(figsize=FIGSIZE)
  tx = datetime.datetime.now()
  e3 = ncB[FIELD].where(((xr.ufuncs.isnan(ncA[FIELD])) & (ncB[FIELD]>=0.0))).values.flatten()
  ne3 = ((xr.ufuncs.isnan(ncA[FIELD])) & (ncB[FIELD]>=0.0)).sum().values
  print(datetime.datetime.now()-tx)
  ax4.set_ylim(bottom=YMIN,top=YMAX)
  ax4.set_title('B >= 0.0 where A was missing\nN = %d\n%f %% of all points' % (ne3,(ne3/nptsanly)*100.0))
  ax4.set_xlabel('Values')
  ax4.set_ylabel('Percent of all values meeting condition')
  e3h = ax4.hist(e3,bins=BINS,density=DENS)
  print(e3h[0])
  print(datetime.datetime.now()-tx)
  fig4.savefig('e3_hg.png')

# **C2** Histogram of zero A values where B is missing
# Q: How many occurrences of A changing from 0.0 to missing data were there?
if C2:
  print("======C2======")
  tx = datetime.datetime.now()
  c2 = ((ncA[FIELD]==0.0) & (xr.ufuncs.isnan(ncB[FIELD]))).sum().values
  print("\n\n\n\n\n")
  print("C2 (A = 0.0, B = missing):")
  #print(datetime.datetime.now()-tx)
  print(c2)
  #print(datetime.datetime.now()-tx)

# **B3** Histogram of missing A values where B is zero
# Q: How many occurrences of A changing from missing data to 0.0 were there?
if B3:
  print("======B3======")
  tx = datetime.datetime.now()
  b3 = ((ncB[FIELD]==0.0) & (xr.ufuncs.isnan(ncA[FIELD]))).sum().values
  print("\n\n\n\n\n")
  print("B3 (A = missing, B = 0.0):")
  #print(datetime.datetime.now()-tx)
  print(b3)
  #print(datetime.datetime.now()-tx)

# **XXB5** Histogram of valid A values where B is zero
# This would essentially show the same thing as B4, except with a bunch of zeros in A.
# B4 would be better to view, because it eliminates 0.0-0.0 values between A and B.

# **XXE2** Histogram of zero A values where B is valid
# This would essentially show the same thing as D2, except with a bunch of zeros in B.
# D2 would be bettwe to view, because it eliminates 0.0-0.0 values between A and B.

# **D2** Histogram of zero A values where B is positive
# Q: If there were 0.0 values in A that became > 0.0 in B, what is the distribution of those B values?
if D2:
  print("======D2======")
  fig5,ax5 = plt.subplots(figsize=FIGSIZE)
  tx = datetime.datetime.now()
  d2 = ncB[FIELD].where(((ncB[FIELD]>0.0) & (ncA[FIELD]==0.0))).values.flatten()
  nd2 = ((ncB[FIELD]>0.0) & (ncA[FIELD]==0.0)).sum().values
  print(datetime.datetime.now()-tx)
  ax5.set_ylim(bottom=YMIN,top=YMAX)
  ax5.set_title('B > 0.0 where A was 0.0\nN = %d\n%f %% of all points' % (nd2,(nd2/nptsanly)*100.0))
  ax5.set_xlabel('Values')
  ax5.set_ylabel('Percent of all values meeting condition')
  d2h = ax5.hist(d2,bins=BINS,density=DENS)
  print(d2h[0])
  print(datetime.datetime.now()-tx)
  fig5.savefig('d2_hg.png')

# **B4** Histogram of positive A values where B is zero
# Q: If there were values > 0.0 in A that became 0.0 in B, what is the distribution of those A values?
if B4:
  print("======B4======")
  fig6,ax6 = plt.subplots(figsize=FIGSIZE)
  tx = datetime.datetime.now()
  b4 = ncA[FIELD].where(((ncA[FIELD]>0.0) & (ncB[FIELD]==0.0))).values.flatten()
  nb4 = ((ncA[FIELD]>0.0) & (ncB[FIELD]==0.0)).sum().values
  print(datetime.datetime.now()-tx)
  ax6.set_ylim(bottom=YMIN,top=YMAX)
  ax6.set_title('A > 0.0 where B is 0.0\nN = %d\n%f %% of all points' % (nb4,(nb4/nptsanly)*100.0))
  ax6.set_xlabel('Values')
  ax6.set_ylabel('Percent of all values meeting condition')
  b4h = ax6.hist(b4,bins=BINS,density=DENS)
  print(b4h[0])
  print(datetime.datetime.now()-tx)
  fig6.savefig('b4_hg.png')

# **C4** Histogram of positive A values where B is missing
# Q: If there were A values that were >0.0 that became missing in B, what is the distribution of those A values?
if C4:
  print("======C4======")
  fig7,ax7 = plt.subplots(figsize=FIGSIZE)
  tx = datetime.datetime.now()
  c4 = ncA[FIELD].where(((ncA[FIELD]>0.0) & (xr.ufuncs.isnan(ncB[FIELD])))).values.flatten()
  nc4 = ((ncA[FIELD]>0.0) & (xr.ufuncs.isnan(ncB[FIELD]))).sum().values
  print(datetime.datetime.now()-tx)
  ax7.set_ylim(bottom=YMIN,top=YMAX)
  ax7.set_title('A > 0.0 where B is missing\nN = %d\n%f %% of all points' % (nc4,(nc4/nptsanly)*100.0))
  ax7.set_xlabel('Values')
  ax7.set_ylabel('Percent of all values meeting condition')
  c4h = plt.hist(c4,bins=BINS,density=DENS)
  print(c4h[0])
  print(datetime.datetime.now()-tx)
  plt.savefig('c4_hg.png')

# **D3** Histogram of missing A, where B is positive
# Q: If there were B values that were >0.0 that were previously missing in A, what is the distribution of those B values?
if D3:
  print("======D3======")
  fig8,ax8 = plt.subplots(figsize=FIGSIZE)
  tx = datetime.datetime.now()
  d3 = ncB[FIELD].where(((ncB[FIELD]>0.0) & (xr.ufuncs.isnan(ncA[FIELD])))).values.flatten()
  nd3 = ((ncB[FIELD]>0.0) & (xr.ufuncs.isnan(ncA[FIELD]))).sum().values
  print(datetime.datetime.now()-tx)
  ax8.set_ylim(bottom=YMIN,top=YMAX)
  ax8.set_title('B > 0.0 where A was missing\nN = %d\n%f %% of all points' % (nd3,(nd3/nptsanly)*100.0))
  ax8.set_xlabel('Values')
  ax8.set_ylabel('Percent of all values meeting condition')
  d3h = plt.hist(d3,bins=BINS,density=DENS)
  print(d3h[0])
  print(datetime.datetime.now()-tx)
  fig8.savefig('d3_hg.png')

# **XXD5** Histogram of valid A where B is positive
# This would show the same thing as D4, except with the addition of a bunch of 0.0 values in A.
# It would be better to use D4 which eliminates 0.0 values in A, or look explicitly at 0.0 values
# in A that became positive which can be seen in D2.

# **XXE4** Histogram of positive A where B is valid
# This would show the same thing as D4, except with the addition of a bunch of 0.0 values in B.
# It would be better to use D4 which eliminates 0.0 values in B, or look explicitly at 0.0 values
# in B that used to be positive in A which can be seen in B4.

# **B2** Count of zero in A and zero in B
# Q: How many occurrences of A and B both equal to 0.0 were there?
if B2:
  print("======B2======")
  tx = datetime.datetime.now()
  b2 = ((ncB[FIELD]==0.0) & (ncA[FIELD]==0.0)).sum().values
  print("\n\n\n\n\n")
  print("B2 (A = 0.0, B = 0.0):")
  print(datetime.datetime.now()-tx)
  print(b2)
  print(datetime.datetime.now()-tx)

# **C3** Count of missing in A and missing in B
# Q: How many occurrences of A and B both equal to missing were there?
if C3:
  print("======C3======")
  tx = datetime.datetime.now()
  c3 = (xr.ufuncs.isnan(ncB[FIELD]) & xr.ufuncs.isnan(ncA[FIELD])).sum().values
  print("\n\n\n\n\n")
  print("C3 (A = missing, B = missing):")
  print(c3)
  print(datetime.datetime.now()-tx)
  print(datetime.datetime.now()-tx)

# Exit here
exit()

# Replace all missing values if requested
#if subMissing:
#  ncA['ICE_PROB'] = ncA[FIELD].fillna(subMissingA)
#  ncB['ICE_PROB'] = ncB[FIELD].fillna(subMissingB)

# Take a difference
#print("\nCOMPUTING DIFFERENCE...")
#if DIFFDIR=='amb':
#  if subMissing:
#    difference = ncA[FIELD].fillna(subMissingA)-ncB[FIELD].fillna(subMissingB)
#  else:
#    difference = ncA[FIELD]-ncB[FIELD]
#else:
#  if subMissing:
#    difference = ncB[FIELD].fillna(subMissingB)-ncA[FIELD].fillna(subMissingA)
#  else:
#    difference = ncB[FIELD]-ncA[FIELD]

# Make sure there are no fill values, there shouldn't be
#missVals = difference.where(xr.ufuncs.isnan(difference),drop=True)

# Clip the difference values
#diffClip = difference.clip(min=CLIPMIN,max=CLIPMAX)
#print(diffClip.max().values)
#print(diffClip.min().values)

# Mask where differences are >= some thresh. Drop any points where they don't meet the thresh.
print("\nMASKING DIFFERENCE...")
#if DIFFTYPE=='all':
#  # Values in dataset A where the differences are what we want
#  diffA = ncA[FIELD].fillna(subMissingA).where(((diffClip>=DIFFTHRESH) | (diffClip<=-1.0*DIFFTHRESH)),drop=True)
#
#  # Values in dataset B where the differences are what we want
#  diffB = ncB[FIELD].fillna(subMissingB).where(((diffClip>=DIFFTHRESH) | (diffClip<=-1.0*DIFFTHRESH)),drop=True)
#  print(diffA.sizes)
#  print(diffB.sizes)
#elif DIFFTYPE=='pos':
#  diffmask = diffClip.where(diffClip>=DIFFTHRESH,drop=True)
#  #diffsave = difference.where(diffClip>=DIFFTHRESH,drop=True)
#elif DIFFTYPE=='neg':
#  diffmask = diffClip.where(diffClip<=-1.0*DIFFTHRESH,drop=True)
#  #diffsave = difference.where(diffClip<=-1.0*DIFFTHRESH,drop=True)
#else:
#  print("\nFATAL! UNKNOWN DIFFTYPE")
#  exit()

# Now that we have diffA and diffB, we can create additional items we want to interrogate:
# 1.  

# Total number of points with any difference
#nPtsDiff = int(np.prod(list(diffmask.sizes.values())))

# Total number of any points
#nPtsVol = int(np.prod(list(diffClip.sizes.values())))

# Set a new figure window
#fig = plt.figure(1, figsize=(22, 15))
fig = plt.figure(1, figsize=(15, 15))

ax1 = plt.subplot(111)
ax1.set_ylim((-1,1))
ax1.set_xlim((-1,1))

plt.scatter(diffA.values.flatten(),diffB.values.flatten(),marker='.',c=['red'])

# TESTS
#plt.hist2d(diffA.values.flatten(),diffB.values.flatten(),bins=[0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0],cmin=1)
#sns.jointplot(data=(diffA.values.flatten(),diffB.values.flatten()),x="A",y="B",kind="hist")

# Plot
#ax1 = plt.subplot(111)
#ax1.set_xticks(ticks=[-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1])

# Histogram
#hcounts, hbins, hbars = xr.plot.hist(diffmask,ax=ax1,bins=[-1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2],rwidth=0.75)

# Compute bin percentages, compared to total points with any differences
#binPercs = (hcounts/nPtsDiff)*100.0

# Add labels to the top of each bar
#binCnt = 0
#for p in hbars:
#  txX = p.get_x()+p.get_width()/2.0
#  txY = p.get_height()+0.05*max(binPercs)
#  ax1.text(txX,txY,'%0.2f' % (binPercs[binCnt]),ha='center')
#  binCnt = binCnt + 1

# Figure title
#fig.suptitle('Differences > %03f %%, N = %15d' % (DIFFTHRESH,int(np.prod(list(diffmask.sizes.values())))))

# Save the image
plt.savefig('diffscatter.png')
