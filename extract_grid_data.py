#!/usr/local/python/bin/python

# Modules
import xarray as xr
import pandas as pd
import sys
import statistics as stat
import numpy as np
from scipy import stats
import subprocess, os, time

############################# User Config #################################

# TODO:
# 1. Spot check in NCL
# 2. Write out min/max values in neighborhood (finalhood) then we can do scoring elsewhere

# Model string
#mstring = "num25_hrrr"
#mid = "hrrr"
mstring = "num4_rap"
mid = "rap"

# Define chunks. Try to keep chunks ~ 1M pts
# Higher numbers = smaller chunks, lower numbers = bigger chunks
#chunks = {'y0':158,'x0':158} # HRRR (40x1059x1799), chunk = (40x158x158)
chunks = {'y0':160,'x0':160} # RAP (39x337x451), chunk = (39x160x160)
#chunks = {'y0':225,'x0':225}

# File with PIREP, matching model, and location info
#matchfile = "/home/dadriaan/projects/sae2019/data/match/hrrr/PIREPShrrr21600posneg.out"
#matchfile = "/home/dadriaan/projects/sae2019/data/match/hrrr/PIREPShrrr10800posneg.out"
#matchfile = "/home/dadriaan/projects/sae2019/data/match/rap/PIREPSrap21600posneg.out"
#matchfile = "/home/dadriaan/projects/sae2019/data/match/rap/PIREPSrap10800posneg.out"
matchfile = "/home/dadriaan/projects/sae2019/data/match/rap/pirepTest.out"
#matchfile = "/home/dadriaan/projects/sae2019/data/match/hrrr/pirepTest.out"

# Variable URL's
probURL = "/var/autofs/mnt/ahmose_d1/dadriaan/projects/sae2019/data/"+mstring+"/data/netcdf/fip/pressure/conus_ICE_PROB"
sevscenURL = "/var/autofs/mnt/ahmose_d1/dadriaan/projects/sae2019/data/"+mstring+"/data/netcdf/fip/diagnostic/conus_SEV_SCENARIO"
spcpURL = "/var/autofs/mnt/ahmose_d1/dadriaan/projects/sae2019/data/"+mstring+"/data/netcdf/fip/diagnostic/conus_SURF_PRECIP"
hgtURL = "/var/autofs/mnt/ahmose_d1/dadriaan/projects/sae2019/data/"+mstring+"/data/netcdf/model/"+mid+"/pressure_derived/conus_HGT"
sevURL = "/var/autofs/mnt/ahmose_d1/dadriaan/projects/sae2019/data/"+mstring+"/data/netcdf/fip/pressure/conus_ICE_SEV"
sldURL = "/var/autofs/mnt/ahmose_d1/dadriaan/projects/sae2019/data/"+mstring+"/data/netcdf/fip/pressure/conus_SLD"
rhURL = "/var/autofs/mnt/ahmose_d1/dadriaan/projects/sae2019/data/"+mstring+"/data/netcdf/model/"+mid+"/pressure_derived/conus_RH"
slwURL = "/var/autofs/mnt/ahmose_d1/dadriaan/projects/sae2019/data/"+mstring+"/data/netcdf/model/"+mid+"/pressure_derived/conus_SLW"
tmpURL = "/var/autofs/mnt/ahmose_d1/dadriaan/projects/sae2019/data/"+mstring+"/data/netcdf/model/"+mid+"/pressure_derived/conus_TMP"
vvelURL = "/var/autofs/mnt/ahmose_d1/dadriaan/projects/sae2019/data/"+mstring+"/data/netcdf/model/"+mid+"/pressure_derived/conus_VVEL"
liqcURL = "/var/autofs/mnt/ahmose_d1/dadriaan/projects/sae2019/data/"+mstring+"/data/netcdf/model/"+mid+"/pressure_derived/conus_LIQ_COND"
icecURL = "/var/autofs/mnt/ahmose_d1/dadriaan/projects/sae2019/data/"+mstring+"/data/netcdf/model/"+mid+"/pressure_derived/conus_ICE_COND"

# Height offset to expand PIREP top/base. If this is zero, then the code works differently.
zOff = 0

# If zOff = 0, use this bubble to search for points around the PIREP within 1000 ft
bubble = 1000

# What's the dz of the dataset (if not constant altitude)?
dz = 25

# Column names in CSV file
colNames = ['id','unixObs','lat','lon','flvl','temp','ibase1','itop1','iint1','ityp1','ibase2','itop2','iint2','ityp2','acft','rawOb','unixFcst','forecast','minI','maxI','minJ','maxJ','homeI','homeJ','corner','npts','file_string']

# Threshold for probability to use for scoring
probthresh = 0.05 # (5%)

# Threshold for probability when scoring severity
probsevthresh = 0.05 # (5%)

# Threshold for scoring SLW
slwthresh = 1.0e-6

# Column names for dataframe of results that will be appended to the dataframe of the input data
vxCols = ['id','file_string','probFYOY','probFNON','probFNOY','probFYON','sevFYOY','sevFNON','sevFNOY','sevFYON','VVELnear','VVELmean','VVELmed','VVELstdev','RHnear','RHmean','RHmed','RHstdev','SLWnear','SLWmean','SLWmed','SLWstdev','TMPnear','TMPmean','TMPmed','TMPstdev','PCPnear','PCPmode','PCPuni','PCPs','SCENnear','SCENmode','SCENuni','SCENs','nPtsHood','nLevHood','nPtsStats','isMulti','levNear','zNear','badPirep','slwFYOY','slwFYON','slwFNOY','slwFNON','probMIN','probMAX','sevMIN','sevMAX','VVELmin','VVELmax','RHmin','RHmax','TMPmin','TMPmax','SLWmin','SLWmax','ICECmin','ICECmax','LIQCmin','LIQCmax','TOTCmin','TOTCmax']

# Debug flag
DEBUG = True

# Boolens for processing
probScore = False
sevScore = False
slwScore = False
vNear = False
vMean = False
vMed = False
vStdev = False

# Output dataframe
outfile = "/home/dadriaan/projects/sae2019/scripts/hrrr10800vx.csv"
#outfile = "/home/dadriaan/projects/sae2019/scripts/hrrr21600vx.csv"
#outfile = "/home/dadriaan/projects/sae2019/scripts/rap10800vx.csv"
#outfile = "/home/dadriaan/projects/sae2019/scripts/rap21600vx.csv"

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

# Find unique forecast times. We can then group PIREPs using this key
unique_fcst = data.unixFcst.unique()

# Group the PIREPS using the forecast time. This creates a new dataframe
groups = data.groupby(data.unixFcst)

# Find unique files also, in case we're curious
unique_files = data.file_string.unique()

# Boolean for initializing output data
init = False
  
# Loop over each group of PIREPs (this will iterate once per file
for name, group in groups:
    
  # Create a dataframe for this group of PIREPs
  vxData = pd.DataFrame(columns=vxCols,index=list(range(0,len(group),1)))
  vxData.astype('object')
  vxcnt = 0
  
  # Print the first file_string (same for entire group)
  if DEBUG:
    print("")
    print("PROCESSING MODEL FILE:")
    print(group.file_string.iloc[0])

  # Open multiple files so we have height info
  ncFiles = ['%s/%s' % (probURL,group.file_string.iloc[0]),'%s/%s' % (hgtURL,group.file_string.iloc[0]),'%s/%s' % (sevURL,group.file_string.iloc[0]),'%s/%s' % (sldURL,group.file_string.iloc[0]),'%s/%s' % (sevscenURL,group.file_string.iloc[0]),'%s/%s' % (spcpURL,group.file_string.iloc[0]),'%s/%s' % (vvelURL,group.file_string.iloc[0]),'%s/%s' % (rhURL,group.file_string.iloc[0]),'%s/%s' % (tmpURL,group.file_string.iloc[0]),'%s/%s' % (slwURL,group.file_string.iloc[0]),'%s/%s' % (icecURL,group.file_string.iloc[0]),'%s/%s' % (liqcURL,group.file_string.iloc[0])]
  ncData = xr.open_mfdataset(ncFiles,chunks=chunks)
  print(ncData)

  # Correct ICE_PROB and ICE_SEV from NaN to 0.0
  ncData['ICE_SEV'] = ncData.ICE_SEV.fillna(0.0)
  ncData['ICE_PROB'] = ncData.ICE_PROB.fillna(0.0)
  ncData['SEV_SCENARIO'] = ncData.SEV_SCENARIO.fillna(-1.0)
  ncData['TOT_COND'] = ncData['ICE_COND']+ncData['LIQ_COND']

  # Need to correct condensate variables?

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
    vxData['file_string'][vxcnt] = group.file_string.iloc[icnt]

    # Print info
    print("")
    print("PROCESSING PIREP:")
    print(group.rawOb.iloc[icnt])
    print("")
    print("PIREP ID:")
    print(group.id.iloc[icnt])

    # Store the PIREP severity for use in Vx
    pSev = group.iint1.iloc[icnt]

    # Map the PIREP sev to FIP sev
    # Group SEV1-SEV2 mixed with SEV2
    # TRC = 1
    # LGT = 2,3
    # MOD = 4,5
    # HVY = 6,7
    # NULL = -1
    # FIPNO = 0
    # FIPTRC = 1
    # FIPLGT = 2
    # FIPMOD = 3
    # FIPHVY = 4
    repNULL = False
    repTRC = False
    repLGT = False
    repMOD = False
    repHVY = False
    if str(pSev) in ['-1']:
      repNULL = True
    if str(pSev) in ['1']:
      repTRC = True
    if str(pSev) in ['2','3']:
      repLGT = True
    if str(pSev) in ['4','5']:
      repMOD = True
    if str(pSev) in ['6','7']:
      repHVY = True

    # Figure out the height bounds for this PIREP
    #pBot = (group.ibase1.iloc[icnt]*100.0)*.3048 # In meters
    #pTop = (group.itop1.iloc[icnt]*100.0)*.3048  # In meters
    pBot = (group.ibase1.iloc[icnt]*100)          # In feet
    pTop = (group.itop1.iloc[icnt]*100)           # In feet
    if DEBUG:
      print("")
      print("PIREP BASE/TOP:")
      print(pBot,pTop)
    if pBot > pTop:
      vxData['badPirep'][vxcnt] = True
      icnt = icnt + 1
      vxcnt = vxcnt + 1
      continue 

    # Subset the dataset at the specific indexes. Returns Dataset object.
    regionDS = ncData.isel(time=[0],y0=slice(group.minJ.iloc[icnt],group.maxJ.iloc[icnt]+1),x0=slice(group.minI.iloc[icnt],group.maxI.iloc[icnt]+1))
   
    # Convert height to feet inside the dataset (modify in place)
    regionDS['HGT'] = regionDS.HGT*3.281 # Convert gpm to feet to match with PIREP
    
    # Find data close to the PIREP in the vertical
    if pBot == pTop:
      isMulti = False
      # Find all the points where the height is within "bubble" distance from PIREP, and drop data outside
      donutBot = regionDS.where((abs(regionDS.HGT-pBot)<bubble),drop=True)
      if donutBot.z0.values.size==0:
        vxData['badPirep'][vxcnt] = True
        icnt = icnt + 1
        vxcnt = vxcnt + 1
        continue
      else:
        finalHood = regionDS.sel(z0=list(range(int(max(donutBot.z0.values)),int(min(donutBot.z0.values)-dz),-dz)))
        pNear = pBot
    else:
      isMulti = True
      # Find all the points where the height is within the "bubble" distance from the top and bottom of the PIREP,
      # and drop data outside
      donutBot = regionDS.where((abs(regionDS.HGT-pBot)<bubble),drop=True)
      donutTop = regionDS.where((abs(regionDS.HGT-pTop)<bubble),drop=True)

      # In the event an SKC PIREP goes to some absurd height, set the top of the neighborhood to the maxZ
      # If neither a top nor a base was found skip it
      if donutTop.z0.values.size==0 and donutBot.z0.values.size==0:
        # Entire PIREP above model range, just skip it
        vxData['badPirep'][vxcnt] = True
        icnt = icnt + 1
        vxcnt = vxcnt + 1
        continue
        #finalHood = regionDS.sel(z0=list(range(int(max(donutBot.z0.values)),int(min(ncData.z0.values)-dz),-dz)))
      # If no top was found but a base was found, use the maxZ as the top
      elif donutTop.z0.values.size==0 and donutBot.z0.values.size>0:   
        # Construct the vertical neighborhood using the pressure levels and dz
        finalHood = regionDS.sel(z0=list(range(int(max(donutBot.z0.values)),int(min(ncData.z0.values)-dz),-dz)))
      # Found a top but no bottom- just skip this unclear why it happens. I think it is if the PIREP base is really
      # high but no model levels are within 1000ft of it, but there are some model levels within 1000ft of the PIREP
      # top.
      elif donutBot.z0.values.size==0 and donutTop.z0.values.size>0:
        vxData['badPirep'][vxcnt] = True
        icnt = icnt + 1
        vxcnt = vxcnt + 1
        continue
        #finalHood = regionDS.sel(z0=list(range(int(max(ncData.z0
      # Both a top and a base were found
      else:
        # Construct the vertical neighborhood using the pressure levels and dz
        finalHood = regionDS.sel(z0=list(range(int(max(donutBot.z0.values)),int(min(donutTop.z0.values)-dz),-dz)))

      # Compute a representative height near the middle of the range between pBot and pTop
      pNear = stat.mean([pBot,pTop]) 
      
    # Find the closest z-level to the PIREP
    # Requires 2 steps: first, use sel to use the actual height values (coordinate indexing) then use isel to
    # get the correct i,j (positional indexing)
    nearSub = ncData.sel(z0=list(range(int(finalHood.z0.values[0]),int(finalHood.z0.values[-1])-dz,-dz)))
    nearestHood = nearSub.isel(time=[0],y0=group.homeJ.iloc[icnt],x0=group.homeI.iloc[icnt])
   
    # Just find the nearest K-level now to the height closest to the pNear value (meanZPirep)
    nearestFinal = nearestHood.isel(z0=[np.argmin(abs((nearestHood.HGT.values*3.281)-pNear))])
    
    # VARIABLES:
    # ncData: full dataset
    # regionDS: subset of full dataset containing all data within the X-Y region around the PIREP at all Z
    # finalHood: subset of the region dataset containing only those Z within 1000ft of PIREP
    # nearestHood: the closest X-Y point to the PIREP with all Z within 1000ft of PIREP
    # nearestFinal: the closest X-Y point to the PIREP and closest Z to the mean of the PIREP heights
 
    # Set some output data
    vxData['isMulti'][vxcnt] = isMulti
    vxData['levNear'][vxcnt] = nearestFinal.z0.values[0]
    vxData['zNear'][vxcnt] = nearestFinal.HGT.values[0][0]*3.281

    # Set some boolean variables to determine if there are any severities matching the PIREP in the neighborhood
    # TRUE = all missing, none of that severity found
    # FALSE = not all missing, at lease one of that severity found
    #
    # applied NOT: TRUE = at least one severity found of this category
    #              FALSE = none of that severity found of this category
    #fipNONE = not xr.ufuncs.isnan(finalHood.ICE_SEV.where(finalHood.ICE_SEV==0.0).values).all()
    fipTRC = not xr.ufuncs.isnan(finalHood.ICE_SEV.where((finalHood.ICE_SEV>0.0) & (finalHood.ICE_SEV<0.25)).values).all()
    fipLGT = not xr.ufuncs.isnan(finalHood.ICE_SEV.where((finalHood.ICE_SEV>=0.25) & (finalHood.ICE_SEV<0.425)).values).all()
    fipMOD = not xr.ufuncs.isnan(finalHood.ICE_SEV.where((finalHood.ICE_SEV>=0.425) & (finalHood.ICE_SEV<0.75)).values).all()
    fipHVY = not xr.ufuncs.isnan(finalHood.ICE_SEV.where((finalHood.ICE_SEV>=0.75) & (finalHood.ICE_SEV<=1.0)).values).all()
    if finalHood.ICE_SEV.values.sum()==0.0:
      fipNULL = True
    else:
      fipNULL = False

    print("")
    print("PIREP/FIP SEV:")
    print(repNULL,repTRC,repLGT,repMOD,repHVY)
    print(fipNULL,fipTRC,fipLGT,fipMOD,fipHVY)

    # SCORING
    # FYOY = Model/FIP Yes, PIREP Yes (Correct Hit) (PODy)
    # FYON = Model/FIP Yes, PIREP No  (False Alarm)
    # FNOY = Model/FIP No, PIREP Yes  (Miss)
    # FNON = Model/FIP No, PIREP No   (Correct Null) (PODn)
    # POFD = 1-PODn "viewed as a limited measure of false alarms, since it is the fraction of the observed "no" events that were incorrectly forecast as "yes".

    # RULES:
    # 1. For probPODy, if any of the ICE_PROB points are above the probpodthresh then it's a hit.
    # 2. For all ICE_PROB > sevprobthresh, see if there's an ICE_SEV that matches the pSEV. If so, then a hit.
    # 3. Keep track of cases where there was a False Alarm (FYON), or where there was a miss (FNOY) for investigation
    # 4. For PODn_ICE_PROB, if the entire volume of ICE_PROB was 0.0 or NaN, and the PIREP is null, it's correct null.
    # 5. For PODn_ICE_SEV, if the entire volume of ICE_SEV was 0.0 or NaN, and PIREP is null, it's correct null.
    # 6. For distributions of MODEL_VAR, use the closest point. Find all non-NaN points in the bubble column. Three answers needed:
    #    a. nearest- just take the MODEL_VAR at this level
    #    XX b. mean- take the mean of all MODEL_VAR at this level within the X-Y neighborhood
    #    XX c. median- take the median of all MODEL_VAR at this level within the X-Y neighborhood
    #    XX d. stdev- take the stdev of all MODEL_VAR at this level within the X-Y neighborhood
    #    e. mean- take the mean of all MODEL_VAR in the entire bubble
    #    f. median- take the median of all MODEL_VAR in the entire bubble
    #    g. stdev- take the stdev of all MODEL_VAR in the entire bubble
    #    QQ: Which is more appropriate? Entire bubble or X-Y neighborhood? --> neighborhood
    #  7. For FIP diagnostic variables- what do we want to do?
    #     a. For SURF_PRECIP, it's likely to vary in the horizontal. Should we just take the NEAREST? What else could we do?
    #        I think we could also take the count of all the various ones in the neighborhood? Most frequent maybe is better,
    #        compared with nearest also.
    #     b. For SEV_SCENARIO, we probably want the same. We probably want the "most frequent" (mode) and then compare that
    #        with NEAREST

    # Append all data to the CSV dataframe!
    if probScore:
      print("PROB")
      # Probability block- 4 values, can only be one. Don't explicitly set the others they will automatically be NaN
      if pSev > 0 and finalHood.ICE_PROB.max() > probthresh:
        # FYOY
        vxData['probFYOY'][vxcnt] = 1
      if pSev < 0 and finalHood.ICE_PROB.max() == 0.0 and finalHood.ICE_PROB.min() == 0.0:    
        # FNON
        vxData['probFNON'][vxcnt] = 1
      if pSev > 0 and finalHood.ICE_PROB.max() < probthresh:
        # FNOY
        vxData['probFNOY'][vxcnt] = 1
      if pSev < 0 and finalHood.ICE_PROB.max() > 0.0:
        # FYON
        vxData['probFYON'][vxcnt] = 1

    if slwScore:
      print("SLW")
      # SLW block- 4 values, can only be one.
      if pSev > 0 and finalHood.SLW.max() > slwthresh:
        # FYOY
        vxData['slwFYOY'][vxcnt] = 1
      if pSev < 0 and finalHood.SLW.max() == 0.0 and finalHood.SLW.min() == 0.0:
        # FNON
        vxData['slwFNON'][vxcnt] = 1
      if pSev > 0 and finalHood.SLW.max() < slwthresh:
        # FNOY
        vxData['slwFNOY'][vxcnt] = 1
      if pSev < 0 and finalHood.SLW.max() > 0.0:
        # FYON
        vxData['slwFYON'][vxcnt] = 1

### --> NEEDS WORK
    # Threshold the ICE_SEV continuous field using FIP thresholds?
    # For FIP, 0.0 = NONE, 0-0.25 = TRC, 0.25-0.425 = LGT, 0.425-0.75 = MOD, 0.75-1.0 = HVY
    # >= previous and <= current
    # i.e., >= 0.0 and <= 0.25, >= 0.25 and <= 0.425, >= 0.425 and <= 0.75, >= 0.75 and <= 1.0
    # -or-, <= 0.0 = NONE, > 0.0 and < 0.25 = TRC, >= 0.25 and < 0.425 = LGT, >= 0.425 and < 0.75 = MOD, >= 0.75 and <= 1.0 = HVY
    
    # The key may be to create 4 separate variables with all sevs masked outside the range for each category but how will this help?

    # Severity block- 4 values, can only be one. Don't explicitly set the others they will automatically be NaN
    # syncSev is the PIREP severity matched to the FIP thresholds.
    # if syncSev is 3, and there's any 3 FIP severity in the neighborhood then it's a hit (FYOY)
    # if syncSev is 3 and there's other severities >0.0 but no 3's, then it's a miss (FNOY)
    # if syncSev is -1 and all severities are 0.0 then it's a correct null (FNON)
    # if syncSev is -1 and there's other severities >0.0, then it's a false alarm (FYON)

    if sevScore:
      print("SEV")
      if (repTRC and fipTRC) or (repLGT and fipLGT) or (repMOD and fipMOD) or (repHVY and fipHVY):
        # FYOY
        vxData['sevFYOY'][vxcnt] = 1
      else:
        if ((repTRC) and (fipLGT or fipMOD or fipHVY and not fipTRC)) or ((repLGT) and (fipTRC or fipMOD or fipHVY and not fipLGT)) or ((repMOD) and (fipLGT or fipTRC or fipHVY and not fipMOD)) or ((repHVY) and (fipLGT or fipTRC or fipMOD and not fipHVY)) or ((repTRC or repLGT or repMOD or repHVY) and (fipNULL)):
          # FNOY
          vxData['sevFNOY'][vxcnt] = 1
      if (repNULL) and (fipNULL):
        # FNON
        vxData['sevFNON'][vxcnt] = 1
      if (repNULL) and not (fipNULL):
        # FYON
        vxData['sevFYON'][vxcnt] = 1
### <-- NEEDS WORK

    # VARnear: val
    if vNear:
      print("NEAR")
      vxData['VVELnear'][vxcnt] = nearestFinal.VVEL.values.flatten()[0]
      vxData['RHnear'][vxcnt] = nearestFinal.RH.values.flatten()[0]
      vxData['SLWnear'][vxcnt] = nearestFinal.SLW.values.flatten()[0]
      vxData['TMPnear'][vxcnt] = nearestFinal.TMP.values.flatten()[0]
    
    # VARmean: val
    if vMean:
      print("MEAN")
      vxData['VVELmean'][vxcnt] = finalHood.VVEL.mean(dim=list(finalHood.VVEL.dims)).values 
      vxData['RHmean'][vxcnt] = finalHood.RH.mean(dim=list(finalHood.RH.dims)).values 
      vxData['SLWmean'][vxcnt] = finalHood.SLW.mean(dim=list(finalHood.SLW.dims)).values 
      vxData['TMPmean'][vxcnt] = finalHood.TMP.mean(dim=list(finalHood.TMP.dims)).values 
    
    # VARmed: val
    if vMed:
      print("MED")
      #vxData['VVELmed'][vxcnt] = finalHood.VVEL.median(dim=list(finalHood.VVEL.dims)).values
      vxData['VVELmed'][vxcnt] = stat.median(finalHood.VVEL.values.flatten())
      vxData['RHmed'][vxcnt] = stat.median(finalHood.RH.values.flatten())
      vxData['SLWmed'][vxcnt] = stat.median(finalHood.SLW.values.flatten())
      vxData['TMPmed'][vxcnt] = stat.median(finalHood.TMP.values.flatten())
    
    # VARstdev: val
    if vStdev:
      print("STDEV")
      vxData['VVELstdev'][vxcnt] = finalHood.VVEL.std(dim=list(finalHood.VVEL.dims)).values
      vxData['RHstdev'][vxcnt] = finalHood.RH.std(dim=list(finalHood.RH.dims)).values
      vxData['SLWstdev'][vxcnt] = finalHood.SLW.std(dim=list(finalHood.SLW.dims)).values
      vxData['TMPstdev'][vxcnt] = finalHood.TMP.std(dim=list(finalHood.TMP.dims)).values

    # Mins/maxes
    print("MIN/MAX")
    vxData['probMIN'][vxcnt] = finalHood.ICE_PROB.min().values
    vxData['probMAX'][vxcnt] = finalHood.ICE_PROB.max().values
    vxData['sevMIN'][vxcnt] = finalHood.ICE_SEV.min().values
    vxData['sevMAX'][vxcnt] = finalHood.ICE_SEV.max().values
    vxData['VVELmin'][vxcnt] = finalHood.VVEL.min().values
    vxData['VVELmax'][vxcnt] = finalHood.VVEL.max().values
    vxData['RHmin'][vxcnt] = finalHood.RH.min().values
    vxData['RHmax'][vxcnt] = finalHood.RH.max().values
    vxData['TMPmin'][vxcnt] = finalHood.TMP.min().values
    vxData['TMPmax'][vxcnt] = finalHood.TMP.max().values
    vxData['SLWmin'][vxcnt] = finalHood.SLW.min().values
    vxData['SLWmax'][vxcnt] = finalHood.SLW.max().values
    vxData['ICECmin'][vxcnt] = finalHood.ICE_COND.min().values
    vxData['ICECmax'][vxcnt] = finalHood.ICE_COND.max().values
    vxData['LIQCmin'][vxcnt] = finalHood.LIQ_COND.min().values
    vxData['LIQCmax'][vxcnt] = finalHood.LIQ_COND.max().values
    vxData['TOTCmin'][vxcnt] = finalHood.TOT_COND.min().values
    vxData['TOTCmax'][vxcnt] = finalHood.TOT_COND.max().values
    
    # PCPnear: val
    print("PCP")
    vxData['PCPnear'][vxcnt] = ncData.SURF_PRECIP.isel(time=[0],z0=[-1],y0=group.homeJ.iloc[icnt],x0=group.homeI.iloc[icnt]).values.flatten()[0]
    
    # PCPuni/PCPmode:
    if min(regionDS.SURF_PRECIP.isel(z0=[-1]).values.flatten()) == max(regionDS.SURF_PRECIP.isel(z0=[-1]).values.flatten()):
      vxData['PCPuni'][vxcnt] = True
      vxData['PCPmode'][vxcnt] = regionDS.SURF_PRECIP.min().round().values
    else:
      vxData['PCPuni'][vxcnt] = False
      #print(np.unique(regionDS.SURF_PRECIP.isel(z0=[-1]).values.flatten().round()))
      vxData['PCPs'][vxcnt] = np.unique(regionDS.SURF_PRECIP.isel(z0=[-1]).values.flatten().round())
      vxData['PCPmode'][vxcnt] = stats.mode(regionDS.SURF_PRECIP.isel(z0=[-1]).values.flatten())[0][0]
    
    # SCENnear: val
    print("SCEN")
    vxData['SCENnear'][vxcnt] = nearestFinal.SEV_SCENARIO.values.flatten().round()[0]
    
    # SCENuni/SCENmode:
    # Since SEV_SCENARIO is initialized to NaN, the entire neighborhood could be NaN. Leave the SCENuni set to NaN if so.
    if not xr.ufuncs.isnan(finalHood.SEV_SCENARIO.values.flatten().round()).all():
      if min(finalHood.SEV_SCENARIO.values.flatten()).round() == max(finalHood.SEV_SCENARIO.values.flatten()).round():
        vxData['SCENuni'][vxcnt] = True
        vxData['SCENmode'][vxcnt] = finalHood.SEV_SCENARIO.min().round().values
      else:
        vxData['SCENuni'][vxcnt] = False
        #print(np.unique(finalHood.SEV_SCENARIO.values.flatten().round()))
        vxData['SCENs'][vxcnt] = np.unique(finalHood.SEV_SCENARIO.values.flatten().round())
        vxData['SCENmode'][vxcnt] = stats.mode(finalHood.SEV_SCENARIO.values.flatten())[0].round()[0]
    
    # nPtsHood: Number of total non-NaN points in the finalHood
    print("NPTS")
    vxData['nPtsHood'][vxcnt] = finalHood.HGT.size
    
    # nLevHood: Number of vertical levels comprising the finalHood
    vxData['nLevHood'][vxcnt] = finalHood.HGT.sizes['z0']
    
    # nPtsStats: Number of valid points in the neighborhood within 1000ft of PIREP
    #            This is total number of points minus all missing (NaN) points
    vxData['nPtsStats'][vxcnt] = vxData['nPtsHood'][vxcnt]-xr.ufuncs.isnan(finalHood.HGT).sum(dim='z0').sum().values

    # Print out a row of the dataframe to see all the data we've collected for this PIREP
    print("")
    print(vxData.loc[[vxcnt]])

    # Exit after one PIREP
    #sys.exit(0)
    
    # When we get here, see if any of the levels are partially NaNs. This could occur where one gridpoint had
    # height values (or a few GP) that were lower than the rest on the level below or above where they all were
    # within 1000 ft
    # Sum in the z-dimension (smush, composite) and then take a sum of the composite to ensure no nan
    chk = xr.ufuncs.isnan(finalHood.HGT.sum(dim='z0',skipna=False).sum(skipna=False)).values
    if chk:
      print("WARNING! FOUND ODDBALL")
      print(finalHood.dims)
      print(finalHood.HGT.values)
      print(finalHood.ICE_PROB.values)
      #sys.exit(1)
    
    # Advance the PIREP counter
    icnt = icnt + 1
    
    # Advance the vxcnt
    vxcnt = vxcnt + 1
  
    # Print time
    print("")
    print("SECONDS TO PROCESS PIREP:")
    print(time.time()-start_time)
    print("")

  # After each model file, write out the dataframe to CSV, appending the input
  if os.path.exists(outfile):
    vxInput = pd.read_csv(outfile)
    vxInput = vxInput.append(vxData,ignore_index=True)
    print(vxInput)
    vxInput.to_csv(outfile,index=False,na_rep="NaN")
  else:
    vxData.to_csv(outfile,index=False,na_rep="NaN") 

  # Close the dataset
  ncData.close()

  # Exit after one model file
  #sys.exit(0)

