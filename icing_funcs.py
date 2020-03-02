#!/usr/local/python3/bin/python
#
# Functions for icing work

import xarray as xr
import os
import time
import pyart
import math
import datetime
import glob
import itertools

from math import pi
from numpy import cos, sin
from scipy.spatial import cKDTree
import numpy as np

# Define a function to load MDV dataset
def load_mdv_dataset(mdvFiles,DEBUG):

  # For each MDV file:
  # 1. Open file
  # 2. Convert to Xarray
  # 3. Merge with existing dataset if existing dataset
  # 4. Delete single variable dataset
  
  st = time.time()
  fcnt = 0
  for f in mdvFiles:
    # Read the current file into an xarray dataset
    tt1 = time.time()
    fileData = pyart.io.read_grid_mdv(f,file_field_names=True,exclude_fields=[]).to_xarray()
    if DEBUG:
      print("READ")
      print(f)
      print(time.time()-tt1)
    if fcnt==0:
      # If it's the first file, just set allData equal to this files' data
      allData = fileData
    else:
      # If it's not the first file, then set allData to a temporary varaible, delete allData, then reset
      # allData to the merge of the previous allData, and the current files' data
      tmpData = allData
      del(allData)
      tt2 = time.time()
      allData = xr.merge([tmpData,fileData])
      if DEBUG:
        print("MERGE")
        print(time.time()-tt2)
      del(tmpData)

    # Print the current allData to see
    #print(allData)

    # Before moving on, delete the current fileData dataset
    fileData.close()
    del(fileData)
    fcnt = fcnt + 1

  # Return the MDV dataset object
  if DEBUG:
    print(time.time()-st)
  allData = allData.rename({'x':'x0','y':'y0','z':'z0'})
  return(allData)

# Load HRRR projection params into a CF-compliant dictionary
def load_hrrr_proj(data):

  cfparams = {}
  cfparams['grid_mapping_name'] = "lambert_conformal_conic"
  cfparams['standard_parallel'] = 38.5
  cfparams['longitude_of_projection_origin'] = -97.5
  cfparams['latitude_of_projection_origin'] = 38.5
  cfparams['earth_radius'] = 6371229.0 # HRRR SPHERE
  cfparams['projection_x_coordinate'] = data.x.values
  cfparams['projection_y_coordinate'] = data.y.values
  #cfparams['ellps'] = 'sphere'
  #cfparams['earth_radius'] = 6378137.0 # EQUATORIAL, OLD MDV
  #cfparams['earth_radius'] = 6370997.0 # PROJ SPHERE
  #cfparams['earth_radius'] = 6371000.0 # MDV SPHERE
  
  return(cfparams)

# Load RAP projection params into a CF-compliant dictionary
def load_rap_proj(data):
  
  cfparams = {}
  cfparams['grid_mapping_name'] = "lambert_conformal_conic"
  cfparams['standard_parallel'] = 25.0
  cfparams['longitude_of_projection_origin'] = -95.0
  cfparams['latitude_of_projection_origin'] = 25.0
  cfparams['earth_radius'] = 6371229.0 # RAP SPHERE
  cfparams['projection_x_coordinate'] = data.x.values
  cfparams['projection_y_coordinate'] = data.y.values
  #cfparams['ellps'] = 'sphere'
  #cfparams['earth_radius'] = 6378137.0 # EQUATORIAL, OLD MDV
  #cfparams['earth_radius'] = 6370997.0 # PROJ SPHERE
  #cfparams['earth_radius'] = 6371000.0 # MDV SPHERE

  return(cfparams)

# Load GOES-16 projection params into a CF_compliant dictionary
def load_goes16_proj(data):

  cfparams = {}
  cfparams['grid_mapping_name'] = "lambert_conformal_conic"
  cfparams['standard_parallel'] = 25.0
  cfparams['longitude_of_projection_origin'] = -95.0
  cfparams['latitude_of_projection_origin'] = 25.0
  cfparams['earth_radius'] = 6371000.0 # MDV SPHERE
  cfparams['projection_x_coordinate'] = data.x.values
  cfparams['projection_y_coordinate'] = data.y.values
  #cfparams['ellps'] = 'sphere'
  #cfparams['earth_radius'] = 6378137.0 # EQUATORIAL, OLD MDV
  #cfparams['earth_radius'] = 6370997.0 # PROJ SPHERE
  #cfparams['earth_radius'] = 6371000.0 # MDV SPHERE

  return(cfparams)

# KDTree for lat/lon searching
class Kdtree_fast(object):
    def __init__(self, lons, lats):
        rad_factor = pi/180.0 # for trignometry, need angles in radians
        self.latvals = lats[:] * rad_factor
        self.lonvals = lons[:] * rad_factor
        self.shape = self.latvals.shape
        clat,clon = cos(self.latvals),cos(self.lonvals)
        slat,slon = sin(self.latvals),sin(self.lonvals)
        clat_clon = clat*clon
        clat_slon = clat*slon
        triples = list(zip(np.ravel(clat*clon), np.ravel(clat*slon), np.ravel(slat)))
        self.kdt = cKDTree(triples)

    def query(self,lat0,lon0):
        rad_factor = pi/180.0
        lat0_rad = lat0 * rad_factor
        lon0_rad = lon0 * rad_factor
        clat0,clon0 = cos(lat0_rad),cos(lon0_rad)
        slat0,slon0 = sin(lat0_rad),sin(lon0_rad)
        dist_sq_min, minindex_1d = self.kdt.query([clat0*clon0,clat0*slon0,slat0])
        iy_min, ix_min = np.unravel_index(minindex_1d, self.shape)
        return iy_min,ix_min

# Distance calculator using tunnel distance
def ll_distance(lat1,lon1,lat2,lon2,erad):
  if erad > 7000.0:
    erad = erad/1000.0
  rad_factor = pi/180.0
  lat1_rad = lat1 * rad_factor
  lon1_rad = lon1 * rad_factor
  lat2_rad = lat2 * rad_factor
  lon2_rad = lon2 * rad_factor
  dX = cos(lat2_rad)*cos(lon2_rad) - cos(lat1_rad)*cos(lon1_rad)
  dY = cos(lat2_rad)*sin(lon2_rad) - cos(lat1_rad)*sin(lon1_rad)
  dZ = sin(lat2_rad) - sin(lat1_rad)
  return(math.sqrt((math.pow(dX,2))+(math.pow(dY,2))+(math.pow(dZ,2)))*erad)

# Return time/date file formats based on a UNIX timestamp
def format_filename(t,fcst,prod,ff):
  dt = datetime.datetime.fromtimestamp(t)
  # Format = YYYYMMDD/g_HHHHHH/f_SSSSSSSS.FF
  if prod=='fip':
    return(dt.strftime('%Y%m%d')+"/g_"+dt.strftime('%H').zfill(2)+"0000/f_"+str(int(fcst)*3600).zfill(8)+"."+ff)
  if prod=='cip':
    # Format = YYYYMMDD/YYYYMMDD_HHMMSS.FF
    if ff=='nc':
      return(dt.strftime('%Y%m%d')+"/"+dt.strftime('%Y%m%d')+"_"+dt.strftime('%H').zfill(2)+dt.strftime('%M').zfill(2)+dt.strftime('%S').zfill(2)+".nc")
    # Format = YYYYMMDD/HHMMSS.FF
    if ff=='mdv':
      return(dt.strftime('%Y%m%d')+"/"+dt.strftime('%H').zfill(2)+dt.strftime('%M').zfill(2)+dt.strftime('%S').zfill(2)+".mdv")
  
# Function to convert a list of obs filenames/paths to unix times  
def obs_file_list_to_unix(flist,url):
  
  # Get a list of unix times for available files
  flTD = [str.split(f,(url+"/"))[-1] for f in flist]

  # Convert to unix time
  return([datetime.datetime.timestamp(datetime.datetime.strptime(f,'%Y%m%d/%H%M%S.mdv')) for f in flTD])

# Function to locate the closest obs file to a time
def find_closest_obs_file(url,ff,obtime,window,dt):

  # Create a window for searching
  tlower = obtime-window
  tupper = obtime+window

  # Create datestrings from tlower and tupper
  tlds = datetime.datetime.fromtimestamp(tlower).strftime('%Y%m%d')
  tuds = datetime.datetime.fromtimestamp(tupper).strftime('%Y%m%d')
  
  # Create an initial file list
  if os.path.exists('%s/%s' % (url,tlds)):
    flist = glob.glob('%s/%s/*.%s' % (url,tlds,ff))
    flist.sort()
  else:
    print("")
    print("WARNING! PATH DOES NOT EXIST IN find_closest_obs_file.")
    print("")
    exit()

  # Find differences
  diffL = [abs(f-obtime) for f in obs_file_list_to_unix(flist,url)]

  # If the tlds!=tuds, then repeat the above but for tuds
  if tlds!=tuds:
    if os.path.exists('%s/%s' % (url,tuds)):
      flist2 = glob.glob('%s/%s/*.%s' % (url,tuds,ff))
      flist2.sort()
    else:
      print("")
      print("WARNING! PATH DOES NOT EXIST IN find_closest_obs_file.")
      print("")
      exit()

    # Find differences
    diffU = [abs(f-obtime) for f in obs_file_list_to_unix(flist2,url)]

    # Return the unix time of the file matching the obtime
    if min(diffU)<min(diffL):
      #return(flist2[diffU.index(min(diffU))])
      if min(diffU)>dt:
        return(-9999.0)
      else:
        return(obs_file_list_to_unix(flist2,url)[diffU.index(min(diffU))])
    else:
      #return(flist[diffL.index(min(diffL))])
      if min(diffL)>dt:
        return(-9999.0)
      else:
        return(obs_file_list_to_unix(flist,url)[diffL.index(min(diffL))])
  else:
    #return(flist[diffL.index(min(diffL))])
    if min(diffL)>dt:
      return(-9999.0)
    else:
      return(obs_file_list_to_unix(flist,url)[diffL.index(min(diffL))])

