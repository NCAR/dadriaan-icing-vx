#!/usr/bin/env python
import xarray as xr
import glob, time, datetime, sys, os, copy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import numpy as np
import matplotlib

from matplotlib.colors import BoundaryNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable

#from AB_funcs import annotate_heatmap
#from AB_funcs import heatmap

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
LEVELS = np.array([0.0,10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0])
YMIN = 0.0
YMAX = 100.0
FIGSIZE = (20,15)
DEBUG = True
FIELD = "POT_SCENARIO"
UNITS = ""

# Severity scenarios (CIP)
#CAT = np.array([1,2,3,4,5,6,7,9,-9])

# Potential scenarios (CIP)
CAT = np.array([0,1,2,3,4,5,6,-9])
potA = ("TS (0)","NP (1)","OS (2)","NBW (3)","BWN (4)","NWNNZ (5)","NWNZ (6)","NA")
potB = potA

# SLD scenarios (CIP)
#CAT = np.array([0,1,2,3,4,5,6,7,-9])

# File paths
# FIRST
#Apath = "/d1/dadriaan/data/versionB/netcdf/conus_ICE_PROB" # BASELINE
#Bpath = "/d1/dadriaan/data/versionC/ral-ifi-cip/main/data_dir/num25_hrrr/data/netcdf/cip/pressure/conus_ICE_PROB" # NOBLEND

# SECOND
#Apath = "/d1/dadriaan/data/versionC/ral-ifi-cip/main/data_dir/num25_hrrr/data/netcdf/cip/pressure/conus_ICE_PROB" # NOBLEND
#Bpath = "/var/autofs/mnt/khaba_d1/dadriaan/data/versionC/ral-ifi-cip/main/data_dir/num25_hrrr/data/netcdf/cip/pressure/conus_ICE_PROB" # RADIUS

# THIRD
Apath = "/var/autofs/mnt/khaba_d1/dadriaan/data/versionC/ral-ifi-cip/main/data_dir/num25_hrrr/data/netcdf/cip/diagnostic/conus_POT_SCENARIO" # RADIUS
Bpath = "/var/autofs/mnt/khaba_d1/dadriaan/data/versionC/ral-ifi-research/glm_test/data_dir/num25_hrrr/data/netcdf/cip/diagnostic/conus_POT_SCENARIO" # GLM

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

# Since we are just aggregating, we can loop and do this file pair by file pair
fcnt = 0
totpts = 0
totcount = np.zeros([len(CAT),len(CAT)])
while fcnt < len(Afiles):

  if DEBUG:
    print("PROCESSING %s" % (os.path.basename(Afiles[fcnt])))

  # Open the files with Xarray
  t1 = datetime.datetime.now()
  print("READING A...")
  #ncA = xr.open_mfdataset(Afiles,combine='by_coords')
  #ncA = xr.open_mfdataset(Afiles,chunks=chunks,combine='by_coords')
  ncA = xr.open_dataset(Afiles[fcnt],chunks=chunks)
  t2 = datetime.datetime.now()
  print(t2-t1)
  print(ncA)
  print("\nREADING B...")
  #ncB = xr.open_mfdataset(Bfiles,combine='by_coords')
  #ncB = xr.open_mfdataset(Bfiles,chunks=chunks,combine='by_coords')
  ncB = xr.open_dataset(Bfiles[fcnt],chunks=chunks)
  t3 = datetime.datetime.now()
  print(t3-t2)
  print(ncB)

  # Get total points in grid for normalization
  ptsingrid = float(ncA.stack(stacked_dim=[...]).dims.get('stacked_dim'))
  print(ptsingrid)

  # Develop an np ndarray to hold the values
  # Y-axis should be "A" and X-axis should be "B"
  # So for each scenario on the Y-axis, we want to say:
  # how many of this scenario were there in
  counts = np.zeros([len(CAT),len(CAT)])

  # Get the number of rows and columns
  nrow, ncol = counts.shape

  # Set an offset
  off = 0.5

  # What about just writing traditional for loops?
  # The code below implies that A is on the y-axis (rows), and B on the x-axis (cols), and is indexed via data[col,row] when looping rows,cols
  # Loop column by column, row by row
  # This is X in A, when Y in B nrows times
  rcnt = 0
  while rcnt < nrow:
    ccnt = 0
    while ccnt < ncol:
      r = CAT[rcnt]
      c = CAT[ccnt]
      #if DEBUG:
      #  print("COMBO A = %d B = %d" % (int(r),int(c)))
      # If we've reached the NA item in the row count, this means we need to check NA A vs positive B
      if r < 0 and c >= 0:
        counts[ccnt,rcnt] = (xr.ufuncs.isnan(ncA.POT_SCENARIO) & (ncB.POT_SCENARIO>=c-off) & (ncB.POT_SCENARIO<=c+off)).sum().values
      # If we've reached the NA item in the col count, this means we need to check NA B vs. positive A
      elif r >= 0 and c < 0:
        counts[ccnt,rcnt] = ((ncA.POT_SCENARIO>=c-off) & (ncA.POT_SCENARIO<=c+off) & (xr.ufuncs.isnan(ncB.POT_SCENARIO))).sum().values
      # If both items are good, then check those cases
      elif r >= 0 and c >= 0:
        counts[ccnt,rcnt] = (((ncA.POT_SCENARIO>=r-off) & (ncA.POT_SCENARIO<=r+off) & (ncB.POT_SCENARIO>=c-off) & (ncB.POT_SCENARIO<=c+off)).sum().values)
      # If both are NA, check that case
      elif CAT[rcnt] < 0 and CAT[ccnt] < 0:
        counts[ccnt,rcnt] = (xr.ufuncs.isnan(ncA.POT_SCENARIO) & xr.ufuncs.isnan(ncB.POT_SCENARIO)).sum().values
      else:
        print("FATAL.")
        exit()

      ccnt = ccnt + 1
    rcnt = rcnt + 1

  # Before moving to the next file pair, increment totcount, totpts, and fcnt
  totpts = totpts + ptsingrid
  totcount = totcount + counts
  fcnt = fcnt + 1

# After we've covered all files, normalize based on total number of points processed
counts = ((totcount/totpts)*100.0)

# Where counts is zero, set to nan to avoid contouring
counts[counts==0.0] = np.nan

# New figure for heatmap
# The code below was adapted from:
# https://matplotlib.org/3.3.3/gallery/images_contours_and_fields/image_annotated_heatmap.html 
fig, ax = plt.subplots(figsize=FIGSIZE)

# Set normalization for colormap
norm = BoundaryNorm(LEVELS,ncolors=plt.get_cmap('YlGn').N,clip=True)

# Set missing values to gray
mycmap = copy.copy(matplotlib.cm.get_cmap("YlGn"))
mycmap.set_bad(color='lightgray')

# Plot the heatmap
data = counts[::-1]
im = ax.imshow(data, norm=norm,  cmap=mycmap)

# Create colorbar axis
#divider = make_axes_locatable(ax)
#cax = divider.append_axes("right", size="5%", pad=0.05)

# Create colorbar
cbar = ax.figure.colorbar(im, ax=ax, ticks=LEVELS,fraction=0.045,pad=0.05)
cbar.ax.set_ylabel('Percent of all points', rotation=-90, va="bottom")

# We want to show all ticks...
row_labels = potA
col_labels = potB
ax.set_xticks(np.arange(data.shape[1]))
ax.set_yticks(np.arange(data.shape[0]))
# ... and label them with the respective list entries.
ax.set_xticklabels(col_labels)
ax.set_yticklabels(row_labels[::-1])

# Let the horizontal axes labeling appear on top.
ax.tick_params(top=False, bottom=True,
               labeltop=False, labelbottom=True)

# Rotate the tick labels and set their alignment.
plt.setp(ax.get_xticklabels(), rotation=0, ha="center",
         rotation_mode="anchor")

# Turn spines off and create white grid.
for edge, spine in ax.spines.items():
    spine.set_visible(False)

ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
ax.tick_params(which="minor", bottom=False, left=False)
ax.set_ylabel('Dataset A')
ax.set_xlabel('Dataset B')

# Get the data used to annotate the map
data = im.get_array()

# Set some keywords for the text labelling
# Set default alignment to center, but allow it to be
# overwritten by textkw.
kw = dict(horizontalalignment="center",
          verticalalignment="center")

# Set the string formatting of the labels
valfmt = matplotlib.ticker.StrMethodFormatter("{x:.4f} %")

# Loop over the data and create a `Text` for each "pixel".
# Change the text's color depending on the data.
texts = []
textcolors=("black", "white")
threshold = im.norm(data.max())/2.
for i in range(data.shape[0]):
  for j in range(data.shape[1]):
    if np.ma.is_masked(data[i,j]):
      kw.update(color='black')
      text = im.axes.text(j, i, 'NA', **kw)
    else:
      #kw.update(color=textcolors[int(im.norm(data[i, j]) > threshold)])
      kw.update(color='black')
      text = im.axes.text(j, i, valfmt(data[i, j], None), **kw)
    texts.append(text)

fig.savefig('catscatter.png')
