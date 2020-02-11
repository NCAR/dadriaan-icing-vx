#!/usr/local/python3/bin/python
#
# Functions for icing work

import xarray as xr
import os
import time
import pyart

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
  return(allData)

