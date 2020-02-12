#!/usr/bin/env python

# Import ConfigMaster
from ConfigMaster import ConfigMaster

# Create a params class
class Params(ConfigMaster):
  defaultParams = """

#!/usr/bin/env python

# DEFAULT PARAMS

####### mstring #######
#
# NAME: mstring
# OPTIONS:
# TYPE: string
# FORMAT:
# DEFAULT:
# DESCRIPTION: Provide a string representing the model and the number of tiles used (e.g. "num25_hrrr", "num4_rap")
#
mstring = "num25_hrrr"

####### m_id #######
#
# NAME: m_id
# OPTIONS:
# TYPE: string
# FORMAT:
# DEFAULT:
# DESCRIPTION: Provide a model ID as a string (e.g. "hrrr", "rap")
#
m_id = "hrrr"

####### murl_pref #######
#
# NAME: url_pref
# OPTIONS:
# TYPE: string
# FORMAT:
# DEFAULT:
# DESCRIPTION: Prefix for each of the URL's used for model data input
#
murl_pref = ""

###### aurl_pref ######
#
# NAME: aurl_pref
# OPTIONS:
# TYPE: string
# FORMAT:
# DEFAULT:
# DESCRIPTION: Prefix for each of the URL's used for algorithm data input
#
aurl_pref = ""

####### ychunks #######
#
# NAME: ychunks
# OPTIONS:
# TYPE: int
# FORMAT:
# DEFAULT:
# DESCRIPTION: Provide the size of each chunk in the y direction
#              HRRR: 158
#              RAP:  160
#
ychunks = 158

####### xchunks #######
#
# NAME: xchunks
# OPTIONS:
# TYPE: int
# FORMAT:
# DEFAULT:
# DESCRIPTION: Provide the size of each chunk in the x direction
#              HRRR: 158
#              RAP:  160
xchunks = 158

####### avars #######
#
# NAME: avars
# OPTIONS:
# TYPE: list
# FORMAT:
# DEFAULT:
# DESCRIPTION: Provide a list of variables we want to read using their XXXX suffix of the conus_XXXX url for
#              algorithm output.
#
avars = ['ICE_PROB','ICE_SEV','SLD','SEV_SCENARIO','POT_SCENARIO','SLD_SCENARIO']

####### mvars #######
#
# NAME: mvars
# OPTIONS:
# TYPE:
# FORMAT:
# DEFAULT:
# DESCRIPTION: Provide a list of variables we want to read using their XXXX suffix of the conus_XXXX url for
#              model output.
#
mvars = ['HGT','RH','SLW','TMP','VVEL','LIQ_COND','ICE_COND']

####### bubble #######
#
# NAME: bubble
# OPTIONS:
# TYPE: int
# FORMAT: units = feet
# DEFAULT: 1000
# DESCRIPTION: Provide the distance to search for model data within the PIREP height (search above/below the PIREP by
#              this amount
#
bubble = 1000

####### dz #######
#
# NAME: dz
# OPTIONS:
# TYPE: int
# FORMAT: units = hPa
# DEFAULT: 25
# DESCRIPTION: delta-Z between vertical levels of the data
#
dz = 25

####### debug #######
#
# NAME: debug
# OPTIONS:
# TYPE: bool
# FORMAT:
# DEFAULT: True
# DESCRIPTION: Boolean for turning on extra debugging info
#
debug = True

####### vNear #######
#
# NAME: vNear
# OPTIONS:
# TYPE: bool
# FORMAT:
# DEFAULT: False
# DESCRIPTION: Pull out the nearest NWP data to the PIREP
#
vNear = False

####### vMean #######
#
# NAME: vMean
# OPTIONS:
# TYPE: bool
# FORMAT:
# DEFAULT: False
# DESCRIPTION: Pull out the mean of the NWP data in the neighborhood
#
vMean = False

####### vMed #######
#
# NAME: vMed
# OPTIONS:
# TYPE: bool
# FORMAT:
# DEFAULT: False
# DESCRIPTION: Pull out the median of the NWP data in the neighborhood
#
vMed = False

####### vStdev #######
#
# NAME: vStdev
# OPTIONS:
# TYPE: bool
# FORMAT:
# DEFAULT: False
# DESCRIPTION: Pull out the standard deviation of the NWP data in the neighborhood
#
vStdev = False

####### vScen #######
#
# NAME: vScen
# OPTIONS:
# TYPE: bool
# FORMAT:
# DEFAULT: False
# DESCRIPTION: Process scenarios or not
#
vScen = False

####### vNWP #######
#
# NAME: vNWP
# OPTIONS:
# TYPE: bool
# FORMAT:
# DEFAULT: False
# DESCRIPTION: Process NWP data or not
#
vNWP = False

####### vAlgo #######
#
# NAME: vAlgo
# OPTIONS:
# TYPE: bool
# FORMAT:
# DEFAULT: False
# DESCRIPTION: Process algorithm data or not
#
vAlgo = False

####### nsecfcst #######
#
# NAME: nsecfcst
# OPTIONS:
# TYPE: int
# FORMAT:
# DEFAULT: 10800
# DESCRIPTION: Number of seconds in the forecast lead time
#
nsecfcst = 10800

####### outfile #######
#
# NAME: outfile
# OPTIONS:
# TYPE: string
# FORMAT:
# DEFAULT: "tmp.csv"
# DESCRIPTION: Full path to the output CSV file. Should be called "tmp.csv"
#
outfile = "tmp.csv"

####### infile #######
#
# NAME: infile
# OPTIONS:
# TYPE: string
# FORMAT:
# DEFAULT: "tmp.in"
# DESCRIPTION: Full path to the input CSV file. Should be called "tmp.in"
#
infile = "tmp.in" 

####### file_format #######
#
# NAME: file_format
# OPTIONS: "mdv" or "netcdf"
# TYPE: string
# FORMAT:
# DEFAULT: "mdv"
# DESCRIPTION: Specify the algorithm/model input file format
#
file_format = "mdv"

"""
