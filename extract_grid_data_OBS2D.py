#!/usr/local/python3/bin/python

# Import params
from extract_params import Params
p = Params()
p.init()

# Modules
import xarray as xr
import pandas as pd
import sys
import statistics as stat
import numpy as np
from scipy import stats
import subprocess, os, time
import datetime

# Add current path
sys.path.append('.')

# Import icing funcs
import icing_funcs as icing

############################# User Config #################################

# Debug flag
DEBUG = p.opt['debug']

# What file format?
# Supported options are MDV or NETCDF
ff = p.opt['file_format']

# Model string
mstring = p.opt['mstring']

# Define chunks. Try to keep chunks ~ 1M pts
# Higher numbers = smaller chunks, lower numbers = bigger chunks
# HRRR_prs (40x1059x1799), chunk = 40x158x158
# RAP_prs (39x337x451), chunk = 39x160x160
if ff=="mdv":
  chunks = {'y0':p.opt['ychunks'],'x0':p.opt['xchunks']}
elif ff=="netcdf":
  chunks = {'y0':p.opt['ychunks'],'x0':p.opt['xchunks']}
else:
  print("")
  print("UNKNOWN FILE FORMAT.")
  sys.exit(1)

# File with PIREP, matching model, and location info
matchfile = p.opt['infile']

# List to hold obs URL's
ourls = []

# Fill up lists
for v in p.opt['ovars']:
  ourls.append(os.path.join(p.opt['ourl_pref'],mstring,'data/%s/satellite/C02' % (ff)))

# Set the variable name
VARNAME = p.opt['ovars'][0]

# Column names in CSV file
colNames = ['id','unixObs','lat','lon','flvl','temp','ibase1','itop1','iint1','ityp1','ibase2','itop2','iint2','ityp2','acft','rawOb','unixFcst','forecast','minI','maxI','minJ','maxJ','homeI','homeJ','corner','npts','mfile_string','file_string','obsfiletime']

# Column names for dataframe of results that will be appended to the dataframe of the input data
vxCols = ['id','file_string','obsfiletime','VARmin','VARmax','VARnear','VARmean','VARmed','VARstdev','nPtsHood','badPirep']

# Output dataframe
outfile = p.opt['outfile']

############################################################################

# Open match file and store as pandas dataframe object
data = pd.read_csv(matchfile,header=None,names=colNames)
#print(data.columns.get_values())
#probPODy = list(range(0,len(data.index),1))
#print(len(probPODy))
#print(len(data.index))
#data['probPODy'] = probPODy
#print(data.columns.get_values())

# Set pandas options
pd.set_option('display.max_columns', 100)

# Create a new dataframe to hold output data
#vxData = pd.DataFrame(columns=vxCols,index=list(range(0,len(data.index),1)))
#vxData.astype('object')

# Read in dataframe from last processed file, and trim the input data to process based on the ID
foundInput = False
if os.path.exists(outfile):
  chkInput = pd.read_csv(outfile)
  #cmd = 'rm -rfv %s' % (outfile)
  #subprocess.call('%s' % (cmd), shell=True, executable='/bin/csh')
  foundInput = True
  #print(chkInput['id'])
  startID = chkInput['id'][len(chkInput)-1]+1
  #print(startID)
  procData = data[startID:].reset_index()
  #print(procData['id'][0])
  del(data)
  data = procData
  #print(data.shape)
  #print(data['id'])
  #sys.exit(0)

# Group the PIREPS using the obsfiletime. This is so we only
# open each obs file one time for all relevant PIREPs
groups = data.groupby(data.obsfiletime)
  
# Loop over each group of PIREPs (this will iterate once per file
for name, group in groups:
    
  # Create a dataframe for this group of PIREPs
  vxData = pd.DataFrame(columns=vxCols,index=list(range(0,len(group),1)))
  vxData.astype('object')
  vxcnt = 0
 
  # Print the first file_string (same for entire group)
  if DEBUG:
    print("")
    print("PROCESSING OBS FILE:")
    print(datetime.datetime.strftime(datetime.datetime.fromtimestamp(group.obsfiletime.iloc[0]),"%Y%m%d %H:%M:%S"))
  vxData['obsfiletime'][vxcnt] = group.obsfiletime.iloc[0]

  # We need to develop the actual file path before loading it
  fname = icing.format_filename(group.obsfiletime.iloc[0],3,"cip",ff)

  # Open multiple files so we have height info
  # List of netCDF files to open
  ncFiles = []
  print("")
  print("LOADING DATA")
  if ff=="mdv":
    for u in ourls:
      ncFiles.append('%s/%s' % (u,fname))
    # Double check files exist
    for f in ncFiles:
      if not os.path.exists('%s' % (f)):
        print("")
        print("FATAL! FILE %s DOES NOT EXIST." % (f))
        print("")
        exit()
    ncData = icing.load_mdv_dataset(ncFiles,DEBUG=True,DROPTIME=False)
    if p.opt['chunk']:
      ncData.chunk(chunks=chunks)
  elif ff=="netcdf":
    for u in ourls:
      ncFiles.append('%s/%s' % (u,fname))
    # Double check files exist
    for f in ncFiles:
      if not os.path.exists('%s' % (f)):
        print("")
        print("FATAL! FILE %s DOES NOT EXIST." % (f))
        print("")
        exit()
    #ncData = xr.open_mfdataset(ncFiles,chunks=chunks,combine='by_coords',compat='override')
    if p.opt['chunk']:
      ncData = xr.open_mfdataset(ncFiles,chunks=chunks,combine='by_coords')
    else:
      ncData = xr.open_mfdataset(ncFiles,combine='by_coords')
  else:
    print("")
    print("UNKNOWN FILE FORMAT.")
    sys.exit(1)

  # Loop over each item in the series
  icnt = 0
  while icnt < len(group.minI):

    # Flag for whether PIREP is good or not
    vxData['badPirep'][vxcnt] = False
  
    # Start time
    start_time = time.time()

    # This is the index into the master pandas DF of the current object in the group being processed
    # We can use this to append info to a new dataframe? Or insert into the original dataframe?
    #print(group.index[icnt])

    # Store input info in the output dataframe
    vxData['id'][vxcnt] = group.id.iloc[icnt]
    #vxData['file_string'][vxcnt] = group.file_string.iloc[icnt]
    vxData['file_string'][vxcnt] = fname

    # Print info
    print("")
    print("PROCESSING PIREP:")
    print((group.rawOb.iloc[icnt]))
    print("")
    print("PIREP ID:")
    print((group.id.iloc[icnt]))

    # Figure out the height bounds for this PIREP
    #pBot = (group.ibase1.iloc[icnt]*100.0)*.3048 # In meters
    #pTop = (group.itop1.iloc[icnt]*100.0)*.3048  # In meters
    pBot = (group.ibase1.iloc[icnt]*100)          # In feet
    pTop = (group.itop1.iloc[icnt]*100)           # In feet
    if DEBUG:
      print("")
      print("PIREP BASE/TOP:")
      print((pBot,pTop))
    if pBot > pTop:
      vxData['badPirep'][vxcnt] = True
      icnt = icnt + 1
      vxcnt = vxcnt + 1
      continue 

    # Subset the dataset at the specific indexes. Returns Dataset object.
    regionDS = ncData.isel(time=[0],y0=slice(group.minJ.iloc[icnt],group.maxJ.iloc[icnt]+1),x0=slice(group.minI.iloc[icnt],group.maxI.iloc[icnt]+1))

    # Get the nearest point to the ob
    nearestDS = ncData.isel(time=[0],y0=group.homeJ.iloc[icnt],x0=group.homeI.iloc[icnt])

    # VARIABLES:
    # ncData: full dataset
    # regionDS: subset of full dataset containing all data within the X-Y region around the PIREP at all Z
    # nearestDS = nearest X-Y point to ob
 
    # Compute statistics/collect data
    vxData['VARnear'][vxcnt] = nearestDS['%s' % (VARNAME)].values.flatten()[0]
    vxData['VARmean'][vxcnt] = regionDS['%s' % (VARNAME)].mean(dim=list(regionDS['%s' % (VARNAME)].dims)).values
    vxData['VARmed'][vxcnt] = stat.median(regionDS['%s' % (VARNAME)].values.flatten())
    vxData['VARstdev'][vxcnt] = regionDS['%s' % (VARNAME)].std(dim=list(regionDS['%s' % (VARNAME)].dims)).values
    vxData['VARmin'][vxcnt] = regionDS['%s' % (VARNAME)].min().values
    vxData['VARmax'][vxcnt] = regionDS['%s' % (VARNAME)].max().values

    # nPtsHood: Number of total non-NaN points in the finalHood
    print("NPTS")
    vxData['nPtsHood'][vxcnt] = regionDS['%s' % (VARNAME)].size
    
    # Print out a row of the dataframe to see all the data we've collected for this PIREP
    print("")
    print((vxData.loc[[vxcnt]]))

    # Exit after one PIREP
    #sys.exit(0)
    
    # Advance the PIREP counter
    icnt = icnt + 1
    
    # Advance the vxcnt
    vxcnt = vxcnt + 1
  
    # Print time
    print("")
    print("SECONDS TO PROCESS PIREP:")
    print((time.time()-start_time))
    print("")

  # After each model file, write out the dataframe to CSV, appending the input
  if os.path.exists(outfile):
    vxInput = pd.read_csv(outfile)
    vxInput = vxInput.append(vxData,ignore_index=True)
    #print(vxInput)
    vxInput.to_csv(outfile,index=False,na_rep="NaN")
  else:
    vxData.to_csv(outfile,index=False,na_rep="NaN") 

  # Close the dataset
  ncData.close()

  # Exit after one model file
  #sys.exit(0)

