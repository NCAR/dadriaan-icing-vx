#
# file: compute_percent_volume.py
#
#
# Description:
# "Over a specified time period, compute the total volume of airspace warned by a particular field given a threshold"
#
# Example: From 01 Feb 2019 - 07 Feb 2019, compute the percent volume airspace warned by FIP ICE_PROB > 0.05
#
# Params: start_time, end_time, URL, field_name, field_thresh, nsecfcst
#
# Vertical matching:
# 1. Need to compute half-heights for volume computation (z-dim)
# 2. Need to factor in topography height on bottom end?
#
import icing_funcs as icing
import xarray as xr
import numpy as np

# Turn off scientific notation
np.set_printoptions(suppress=True)

################# CONF
# Set the number of seconds in the forecast lead
nsecfcst = 10800

# Set the dZ value of the vertical dimension
dZ = -25.0 # PRESSURE
#dZ = 500.0 # FLIGHT
#dZ = 1.0 # NATIVE
xv = 1100
yv = 500

# Variable
var = ['ICE_PROB']

# Threshold
thresh = 0.55

# Input file
f = "/d1/dadriaan/projects/sae2019/data/num25_hrrr/data/mdv/fip/pressure/conus_ICE_PROB/20190201/g_000000/f_00010800.mdv"
fz = "/d1/dadriaan/projects/sae2019/data/num25_hrrr/data/mdv/model/hrrr/pressure_derived/conus_HGT/20190201/g_000000/f_00010800.mdv"
#fz = "/d1/dadriaan/projects/sae2019/data/num25_hrrr/data/mdv/model/hrrr/native/conus_HGT/20190201/g_000000/f_00010800.mdv"

# Flag for whether to include algorithm output or not
Algo = True

###################

# Read data
if Algo:
  data = icing.load_mdv_dataset([f,fz],True,False)
else:
  data = icing.load_mdv_dataset([f],True,False)

if Algo:
  # Correct NA data in certain variables
  data['ICE_PROB'] = data.ICE_PROB.fillna(0.0)

  # Correct any INT16 values that are <0 to 0.0, in fields that shouldn't have negative data
  data['ICE_PROB'] = data.ICE_PROB.where(data.ICE_PROB>0.0,0.0)

  # Back out the icing potential value
  data['ICE_POT'] = data['ICE_PROB']/((-0.033*(nsecfcst/3600)+0.84))

# Create a DX and DY array, with the edges set to dx/2 and dy/2
data['DX'] = xr.zeros_like(data['ICE_PROB'])
data['DY'] = xr.zeros_like(data['ICE_PROB'])
data['DX'] = xr.where(data['DX']<1.0,3.0,0.0)
data['DY'] = xr.where(data['DY']<1.0,3.0,0.0)
data['DX'].loc[dict(x0=data.coords["x0"].isel(x0=0))] = data['DX'].loc[dict(x0=data.coords["x0"].isel(x0=0))]/2.0
data['DX'].loc[dict(x0=data.coords["x0"].isel(x0=data.dims['x0']-1))] = data['DX'].loc[dict(x0=data.coords["x0"].isel(x0=data.dims['x0']-1))]/2.0
data['DY'].loc[dict(y0=data.coords["y0"].isel(y0=0))] = data['DY'].loc[dict(y0=data.coords["y0"].isel(y0=0))]/2.0
data['DY'].loc[dict(y0=data.coords["y0"].isel(y0=data.dims['y0']-1))] = data['DY'].loc[dict(y0=data.coords["y0"].isel(y0=data.dims['y0']-1))]/2.0

# Compute the heights of each grid cell to use in the VOL calculation
data['DZ'] = data['HGT'].differentiate('z0',1,)*dZ

# Correct the lowest value and the uppermost value
data['DZ'].loc[dict(z0=data['z0'][0])] = data['DZ'].loc[dict(z0=data['z0'][0])]/2.0
data['DZ'].loc[dict(z0=data['z0'][data.dims['z0']-1])] = data['DZ'].loc[dict(z0=data['z0'][data.dims['z0']-1])]/2.0

# Define a custom function for computing the exponent
def xrpower(a,b):
  func = lambda x, y: np.power(x,y)
  return(xr.apply_ufunc(func,a,b))

# Compute the VOL
# This is dZ*dX*dY
data['VOL'] = (data['DZ']/1000.0*data['DX']*data['DY'])

# DEBUGGING
#print(data['HGT'].isel(x0=xv,y0=yv))
#print(data['DZ'].isel(x0=xv,y0=yv))
#print(data['VOL'].isel(x0=xv,y0=yv))
#print(data['ICE_POT'].isel(x0=xv,y0=yv))
#print(data['VOL'].isel(x0=xv,y0=yv).where(data['ICE_POT'].isel(x0=xv,y0=yv)>=thresh).sum('z0'))

# Print the total VOL sum over entire grid for var and thresh
print("")
print("TOTAL VOLUME WARNED FOR %s of >= %f in km^3: " % (var,thresh))
print(data['VOL'].where(data[var]>=thresh).sum(data.dims))
