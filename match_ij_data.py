#!/usr/local/python3/bin/python
#
# File: match_ij_data.py
#
# Author: D. Adriaansen
#
# Date: 14 Feb 2020
#
# Purpose: match observations to grid cells
#
# Notes:
#__________________________________________________________________

# Import params
from match_ij_params import Params
p = Params()
p.init()

# Load Python libraries
import sys, os
import pandas as pd
import pyart
import metpy
import time
import calendar
import matplotlib.pyplot as plt
import datetime

# Add current path
sys.path.append('.')

# Import icing funcs
import icing_funcs as icing

DEBUG = True

start_time = calendar.timegm(time.gmtime())

################################### User Config ###################################

# MDV Grid
mdv = p.opt['mdv_file']

# Point data file
ptobs = p.opt['point_obs_file']

# Output file path
out = p.opt['outfile']

# Number of points per side on the square of your neighborhood
# For HRRR, 8 = 8x8 = 64 pts = 24 km^2 area
# For RAP, 2 = 2x2 = 4 pts = 26 km^2 area
nside = p.opt['nside']

# Offset, to start at the largest bounding square
offset = (nside/2)-1

# How many total points in the neighborhood
npts = nside*nside

# What forecast lead time in seconds do we want to use for matching model data?
nsecfcst = p.opt['fcst_lead_seconds']

# What time delta do we want for matching model data to obs (seconds)?
matchdt = p.opt['match_offset_seconds']

# Do we want to match before, after, or both from the valid time?
matchwindow = p.opt['match_window_method']

# Number of minutes to match valid time
# For FIP, set to 0. For CIP, use minutes after the hour that CIP was run
matchmins = p.opt['match_time_minutes']

# Name of model
model_name = p.opt['model_name']

#######################################################################################

print("")
print("WARNING! This code expects a model grid with the lower left corner as point (0,0)")
print("")

# Clean up if output file exists
try:
  os.remove(out)
except OSError:
  pass

# Get the column names
column_names = p.opt['input_csv_columns']

# Open input file
df = pd.read_csv(ptobs,header=None,names=column_names,index_col=False,skiprows=p.opt['skipNrows'])

# Rename lat/lon if needed, and also adjust the column names to match for writing output
if p.opt['latlon_colmap'] != {}:
  df = df.rename(columns=p.opt['latlon_colmap'])
  for k in list(p.opt['latlon_colmap'].keys()):
    column_names.remove(k)
  for v in list(p.opt['latlon_colmap'].values()):
    column_names.append(v)

# Drop anywhere lat/lon are missing
df = df.loc[df[df['lat']>-9000.0].index]

# Handle time. If unix_time is provided, then use it. Otherwise construct a unix_time column using the dtinfo conf option
if 'unix_time' not in p.opt['input_csv_columns']:
  print("\nNO UNIX TIME FOUND. CONSTRUCTING USING dtinfo PARAM.")

  # Get the list of column names containing the date info
  keylist = list(p.opt['dtinfo'].keys())

  # Get the column formatting for each date info column
  valuelist = list(p.opt['dtinfo'].values())

  # Create the UNIX_TIME by doing:
  # 1. Apply a lambda function to join each column with date info together
  # 2. Convert it to a datetime object using the formatting specified for each column by the user
  # 3. Subtract off the EPOCH time (1970-1-1 00:00:00) and then floor divide by 1s per the Pandas documentation here:
  #    https://pandas.pydata.org/pandas-docs/stable/user_guide/timeseries.html#from-timestamps-to-epoch
  df['unix_time'] = (pd.to_datetime(df[keylist].apply(lambda x: ''.join(x), axis=1),format=''.join(valuelist)) - pd.Timestamp('1970-1-1')) // pd.Timedelta('1s')

  # Add the new column to the column names
  column_names.append('unix_time')

# Open the grid file
grid = pyart.io.read_grid_mdv(mdv,file_field_names=True,exclude_fields=[])
gds = grid.to_xarray()
del(grid)

# Return the projection information from the icing functions
if model_name=='rap':
  proj_params = icing.load_rap_proj(gds)
elif model_name=='hrrr':  
  proj_params = icing.load_hrrr_proj(gds)
elif model_name=='goes16':
  proj_params = icing.load_goes16_proj(gds)
else:
  print("")
  print("MODEL NOT SUPPORTED")
  print("")
  exit()

# Assign the coordinate reference system using MetPy
gds = gds.metpy.assign_crs(cf_attributes=proj_params)

# Assign the 2D latitude/longitude to the grid
gds = gds.metpy.assign_latitude_longitude(force=True)
grid_lat = gds['latitude'].values
grid_lon = gds['longitude'].values

# Set up the KDTree
ll = icing.Kdtree_fast(gds['longitude'].values,gds['latitude'].values)

# Get sizes of dimensions
numX = gds.dims['x']
numY = gds.dims['y']

# Find the ij's of all the PIREPs
# This will return a list of values like [(i1,j1),(i2,j2)]
# Using list comprehension
#ijs = [ll.query(lat,lon)[0] for lat, lon in zip(df['lat'],df['lon'])]
#print(ijs)
#exit()

# Create a timestamp from the user supplied max_obs_time and filter the obs based on that
obs_time_max = datetime.datetime.timestamp(datetime.datetime.strptime(p.opt['max_obs_time'],"%Y-%m-%d %H:%M:%S"))
obs_time_min = datetime.datetime.timestamp(datetime.datetime.strptime(p.opt['min_obs_time'],"%Y-%m-%d %H:%M:%S"))

# Compute the date attributes of the PIREP time
# Using list comprehension
df['year'] = [int(time.gmtime(x)[0]) for x in df['unix_time']]
df['month'] = [int(time.gmtime(x)[1]) for x in df['unix_time']]
df['day'] = [int(time.gmtime(x)[2]) for x in df['unix_time']]
df['hour'] = [int(time.gmtime(x)[3]) for x in df['unix_time']]
df['mins'] = [int(time.gmtime(x)[4]) for x in df['unix_time']]

# Save a histogram of PIREP minutes
#fig,ax = plt.subplots()
#df.hist('mins',ax=ax,bins=24)
#fig.savefig('mins.png')

# Compute the tlower. This is the top of the beginning of the hour containing the PIREP
# Example: PIREP time = 1137, tlower = 1100
#
# Use list comprehension
df['tlower'] = [calendar.timegm((y,m,d,h,matchmins,0)) for y, m, d, h in zip(df['year'],df['month'],df['day'],df['hour'])]

# Create a mask to find any times where tlower is > PIREP time. This will only occur if matchmins is > 0 (mainly for CIP)
mask1 = df['tlower']>df['unix_time']
if len(df[mask1])>0:
  print("")
  print("CHANGING TLOWER FOR %02d PIREPS" % (int(len(df[mask1]))))
  print("")
  df.loc[df[mask1].index,('tlower')] = df.loc[df[mask1].index,('tlower')]-3600

# Compute the tupper. This is the top of the next hour after the hour contining the PIREP
# Example: PIREP time = 1137, tupper = 1200
#
# Use tlower and add one hour
df['tupper'] = df['tlower']+3600

# Create a timedeltas to use when matching
# This is the distance between the PIREP time and the lower time bound. It should always be positive.
df['timedelta_low'] = df['unix_time']-df['tlower']
df['timedelta_upp'] = df['tupper']-df['unix_time']

# Initialize the match time to the tlower
df['tmatch'] = df['tlower']

# Matching methods
# 1. "both" - create a mask finding all PIREPs between matchdt mins before and matchdt mins after top of hour
# 2. "after" - create a mask finding all PIREPs between top of hour and matchdt mins after top of hour
# 3. "before" - create a mask finding all PIREPs between top of hour and matchdt mins before top of hour
if matchwindow=="both":
  print("")
  print("MATCHING BOTH")
  # Since we initialized to tlower, only reset the tmatch where it's >= matchmins, which means it matches to the next time
  mask2 = df['timedelta_low']>=matchdt
  df.loc[df[mask2].index,('tmatch')] = df.loc[df[mask2].index,('tupper')]
  #fdf = df
  #groups = df.groupby(df['tmatch']>0)
if matchwindow=="after":
  print("")
  print("MATCHING AFTER")
  # Only save PIREPs between tlower and tlower + matcthdt mins
  #mask3 = df['timedelta_low']<matchdt
  #fdf = df[mask3]
  #groups = df.groupby(df['timedelta_low']<matchdt)
  df = df[df['timedelta_low']<matchdt]
if matchwindow=="before":
  print("")
  print("MATCHING BEFORE")
  # Only save PIREPs between tupper and tupper - matchdt mins
  #mask4 = df['timedelta_upp']<matchdt
  #fdf = df[mask4]
  #groups = df.groupby(df['timedelta_upp']<matchdt)
  df = df[df['timedelta_upp']<matchdt]
  df.loc[df.index,('tmatch')] = df['tupper']

if DEBUG:
  print(df[['unix_time','tlower','tupper','timedelta_low','timedelta_upp','tmatch']])
  print(len(df))

# Filter the PIREPs based on the tmatch time. This is so we don't get/use model/obs data outside the time window
df = df.loc[df[df['tmatch']>=obs_time_min].index]
df = df.loc[df[df['tmatch']<=obs_time_max].index]

# Set the filename stuff
df.loc[df.index,('mfile_string')] = [icing.format_filename((t-nsecfcst),nsecfcst/3600,'fip','nc') for t in df['tmatch']]
df.loc[df.index,('file_string')] = [icing.format_filename(t,nsecfcst/3600,'cip','nc') for t in df['tmatch']]

# Set the obsfiletime
if model_name in ['goes16']:
  df.loc[df.index,('obsfiletime')] = [icing.find_closest_obs_file(p.opt['data_url'],'mdv',t,1800,900) for t in df['unix_time']]
  # Filter out any PIREPs that weren't matched to an obs file (obs dataset file doesn't exist)
  df = df.loc[df[df['obsfiletime']>0].index]
else:
  df.loc[df.index,('obsfiletime')] = 0.0
#df.loc[df.index,('obsfiletime')] = 0.0

# Use list comprehension and add columns for i,j of PIREP
ijs = [ll.query(lat,lon) for lat, lon in zip(df['lat'],df['lon'])]
df.loc[df.index,('home_i')] = list(zip(*ijs))[1]
df.loc[df.index,('home_j')] = list(zip(*ijs))[0]
# Make sure i,j are integers
df = df.astype({'home_i':'int32','home_j':'int32'})
df.loc[df.index,('home_la')] = [grid_lat[j][i] for j, i in zip(df['home_j'],df['home_i'])]
df.loc[df.index,('home_lo')] = [grid_lon[j][i] for j, i in zip(df['home_j'],df['home_i'])]

if DEBUG:
  print(df[['home_la','home_lo']])

# Add eight new columns, for UR, LR, UL, LL i,j respectively
df.loc[df.index,('UR_i')] = df['home_i']+1
df.loc[df.index,('UR_j')] = df['home_j']+1
df.loc[df.index,('LR_i')] = df['home_i']+1
df.loc[df.index,('LR_j')] = df['home_j']-1
df.loc[df.index,('UL_i')] = df['home_i']-1
df.loc[df.index,('UL_j')] = df['home_j']+1
df.loc[df.index,('LL_i')] = df['home_i']-1
df.loc[df.index,('LL_j')] = df['home_j']-1

# Filter out the dataframe again, by whether any of the corners moved off the grid
df = df[(df['UR_i']>=0) & (df['UR_i']<numX) & (df['UR_j']>=0) & (df['UR_j']<numY)]
df = df[(df['LR_i']>=0) & (df['LR_i']<numX) & (df['LR_j']>=0) & (df['LR_j']<numY)]
df = df[(df['UL_i']>=0) & (df['UL_i']<numX) & (df['UL_j']>=0) & (df['UL_j']<numY)]
df = df[(df['LL_i']>=0) & (df['LL_i']<numX) & (df['LL_j']>=0) & (df['LL_j']<numY)]

# Add the latitude/longitude of the corner points
df.loc[df.index,('UR_la')] = [grid_lat[j][i] for j, i in zip(df['UR_j'],df['UR_i'])]
df.loc[df.index,('UR_lo')] = [grid_lon[j][i] for j, i in zip(df['UR_j'],df['UR_i'])]
df.loc[df.index,('LR_la')] = [grid_lat[j][i] for j, i in zip(df['LR_j'],df['LR_i'])]
df.loc[df.index,('LR_lo')] = [grid_lon[j][i] for j, i in zip(df['LR_j'],df['LR_i'])]
df.loc[df.index,('UL_la')] = [grid_lat[j][i] for j, i in zip(df['UL_j'],df['UL_i'])]
df.loc[df.index,('UL_lo')] = [grid_lon[j][i] for j, i in zip(df['UL_j'],df['UL_i'])]
df.loc[df.index,('LL_la')] = [grid_lat[j][i] for j, i in zip(df['LL_j'],df['LL_i'])]
df.loc[df.index,('LL_lo')] = [grid_lon[j][i] for j, i in zip(df['LL_j'],df['LL_i'])]

# Compute the distance between the PIREP lat/lon and the latitude/longitude of the corners
df.loc[df.index,('UR_dist')] = [icing.ll_distance(hla,hlo,cla,clo,proj_params['earth_radius']) for hla, hlo, cla, clo in zip(df['lat'],df['lon'],df['UR_la'],df['UR_lo'])]
df.loc[df.index,('LR_dist')] = [icing.ll_distance(hla,hlo,cla,clo,proj_params['earth_radius']) for hla, hlo, cla, clo in zip(df['lat'],df['lon'],df['LR_la'],df['LR_lo'])]
df.loc[df.index,('UL_dist')] = [icing.ll_distance(hla,hlo,cla,clo,proj_params['earth_radius']) for hla, hlo, cla, clo in zip(df['lat'],df['lon'],df['UL_la'],df['UL_lo'])]
df.loc[df.index,('LL_dist')] = [icing.ll_distance(hla,hlo,cla,clo,proj_params['earth_radius']) for hla, hlo, cla, clo in zip(df['lat'],df['lon'],df['LL_la'],df['LL_lo'])]

# Compute a minimum distance of X
# If UR_dist < LR_dist & UL_dist & LL_dist then ENCAPS = "UR", corner = "LL"
# If LR_dist < UL_dist & LL_dist & UR_dist then ENCAPS = "LR", corner = "UL"
# If UL_dist < LL_dist & UR_dist & LR_dist then ENCAPS = "UL", corner = "LR"
# if LL_dist < UR_dist & LR_dist & UL_dist then ENCAPS = "LL", corner = "UR"
df.loc[df[(df['UR_dist']<df['LR_dist'])&(df['UR_dist']<df['UL_dist'])&(df['UR_dist']<df['LL_dist'])].index,('corner')] = "LL"
df.loc[df[(df['LR_dist']<df['UL_dist'])&(df['LR_dist']<df['LL_dist'])&(df['LR_dist']<df['UR_dist'])].index,('corner')] = "UL"
df.loc[df[(df['UL_dist']<df['LL_dist'])&(df['UL_dist']<df['UR_dist'])&(df['UL_dist']<df['LR_dist'])].index,('corner')] = "LR"
df.loc[df[(df['LL_dist']<df['UR_dist'])&(df['LL_dist']<df['LR_dist'])&(df['LL_dist']<df['UL_dist'])].index,('corner')] = "UR"

# Compute the minI, maxI, minJ, maxJ for the bounding box
df.loc[df[df['corner']=="LR"].index,('minI')] = (df['home_i']-1)-offset # Left boundary
df.loc[df[df['corner']=="LR"].index,('minJ')] = (df['home_j']-offset)   # Bottom boundary
df.loc[df[df['corner']=="LR"].index,('maxI')] = (df['minI']+nside)-1    # Right boundary
df.loc[df[df['corner']=="LR"].index,('maxJ')] = (df['minJ']+nside)-1    # Top boundary

df.loc[df[df['corner']=="UR"].index,('minI')] = (df['home_i']-1)-offset # Left boundary
df.loc[df[df['corner']=="UR"].index,('minJ')] = (df['home_j']-1)-offset # Bottom boundary
df.loc[df[df['corner']=="UR"].index,('maxI')] = (df['minI']+nside)-1    # Right boundary
df.loc[df[df['corner']=="UR"].index,('maxJ')] = (df['minJ']+nside)-1    # Top boundary

df.loc[df[df['corner']=="LL"].index,('minI')] = (df['home_i']-offset)   # Left boundary
df.loc[df[df['corner']=="LL"].index,('minJ')] = (df['home_j']-1)-offset # Bottom boundary
df.loc[df[df['corner']=="LL"].index,('maxI')] = (df['minI']+nside)-1    # Right boundary
df.loc[df[df['corner']=="LL"].index,('maxJ')] = (df['minJ']+nside)-1    # Top boundary

df.loc[df[df['corner']=="UL"].index,('minI')] = (df['home_i']-offset) # Left boundary
df.loc[df[df['corner']=="UL"].index,('minJ')] = (df['home_j']-offset)   # Bottom boundary
df.loc[df[df['corner']=="UL"].index,('maxI')] = (df['minI']+nside)-1    # Right boundary
df.loc[df[df['corner']=="UL"].index,('maxJ')] = (df['minJ']+nside)-1    # Top boundary

# Convert i,j to integer
df = df.astype({'minI':'int32','maxI':'int32','minJ':'int32','maxJ':'int32'})

# Ensure everything is inside the grid
df = df[(df['minI']>=0) & (df['minI']<numX) & (df['minJ']>=0) & (df['minJ']<numY)]
df = df[(df['maxI']>=0) & (df['maxI']<numX) & (df['maxJ']>=0) & (df['maxJ']<numY)]

# Compute some statistics about the neighborhood
df.loc[df.index,('totPts')] = ((df['maxI']-df['minI'])+1)*((df['maxJ']-df['minJ'])+1)

# Add some missing stuff here
#df['mfile_string'] = ""
#df['file_string'] = ""
df['nsecfcst'] = nsecfcst

# The columns we want by default from the matching code are:
# ['tmatch','nsecfcst','minI','maxI','minJ','maxJ','home_i','home_j','corner','totPts','mfile_string','file_string','obsfiletime']

# Prepend the user "output_csv_columns" to this list
default_out_cols = ['tmatch','nsecfcst','minI','maxI','minJ','maxJ','home_i','home_j','corner','totPts','mfile_string','file_string','obsfiletime']
all_out_cols = p.opt['output_csv_columns']+default_out_cols

# Reset the index
df = df.reset_index(drop=True)

# Write out
df.to_csv(out,header=False,columns=all_out_cols)
print("")
print("PROCESSING TOOK: "+str(calendar.timegm(time.gmtime())-start_time)+" SECONDS.")
