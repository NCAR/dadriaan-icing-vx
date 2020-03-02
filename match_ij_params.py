#!/usr/bin/env python

# Import ConfigMaster
from ConfigMaster import ConfigMaster

# Create a params class
class Params(ConfigMaster):
  defaultParams = """

#!/usr/bin/env python

# DEFAULT PARAMS

####### mdv_file #######
#
# NAME: mdv_file
# OPTIONS:
# TYPE: string
# FORMAT:
# DEFAULT:
# DESCRIPTION: Provide the path to an MDV file containing the grid information to use for matching to the observations
#
mdv_file = ""

####### point_obs_file #######
#
# NAME: point_obs_file
# OPTIONS:
# TYPE: string
# FORMAT:
# DEFAULT:
# DESCRIPTION: Provide the full path to a CSV file containing point observations. These
#              files are expected to have a certain format, which is the "decoded_pirep" format
#              but comma-delimited not space-delimited. These types of files can be created
#              by running parsePirep.pl, PirepCsv2Spdb, then pirep_spdb2decoded.py.
#
point_obs_file = ""

####### outfile #######
#
# NAME: outfile
# OPTIONS:
# TYPE: string
# FORMAT:
# DEFAULT:
# DESCRIPTION: Path to output file. CSV format.
#
outfile = ""

###### nside ######
#
# NAME: nside
# OPTIONS:
# TYPE: int
# FORMAT:
# DEFAULT: 2
# DESCRIPTION: Number of points per side on the square of the neighborhood.
#              Example: nside = 10, square = 10x10 number of points at grid resolution
#
nside = 2

####### fcst_lead_seconds #######
#
# NAME: fcst_lead_seconds
# OPTIONS:
# TYPE: int
# FORMAT:
# DEFAULT: 10800
# DESCRIPTION: Provide the forecast lead time in seconds of the data being used
#
fcst_lead_seconds = 10800

####### match_offset_seconds #######
#
# NAME: match_offset_seconds
# OPTIONS:
# TYPE: int
# FORMAT:
# DEFAULT: 1800
# DESCRIPTION: Provide the offset in seconds from the match time (valid time)
#
match_offset_seconds = 1800

####### match_window_method #######
#
# NAME: match_window_method
# OPTIONS:
# TYPE: string
# FORMAT:
# DEFAULT: 
# DESCRIPTION: Provide the matching window method. This will use match_offset_seconds from the match time (valid time)
#              and either look before, after, or both. If both, match_offset_seconds is used equally before and
#              after. Example: match_offset_seconds = 1800, match_window_method = "both",
#              match window = [valid_time-1800:valid_time+1800]
#              Typically, for FIP we use "both", and for CIP we use "after"
#
match_window_method = "before"

####### match_time_minutes #######
#
# NAME: match_time_minutes
# OPTIONS:
# TYPE:
# FORMAT:
# DEFAULT: 0
# DESCRIPTION: Provide the minutes to be used for the match time, if not the top of the hour.
#              If top of the hour, set to 0.
#              For FIP, set to 0. For CIP, set to number of minutes after the top of the hour CIP was run.
#
match_time_minutes = 0

####### input_csv_columns #######
#
# NAME: input_csv_columns
# OPTIONS:
# TYPE: list
# FORMAT:
# DEFAULT: []
# DESCRIPTION: Provide a list of column names for the input CSV file.
#              
input_csv_columns = ['unix_time','lat','lon','flvl','cbase1','cvg1','ctop1','cbase2','cvg2','ctop2','vis','obx','temp','wdir','wspd','ibase1','itop1','iint1','ityp1','ibase2','itop2','iint2','ityp2','tbase1','ttop1','tint1','ttyp1','tfreq1','tbase2','ttop2','tint2','ttyp2','tfreq2','actype','rawrep']

####### model_name #######
#
# NAME: model_name
# OPTIONS: "rap" or "hrrr"
# TYPE: string
# FORMAT:
# DEFAULT: "hrrr"
# DESCRIPTION: Provide the model name as a string 
#
model_name = "hrrr"

####### max_obs_time #######
#
# NAME: max_obs_time
# OPTIONS:
# TYPE: string
# FORMAT: YYYY-MM-DD HH:MM:SS
# DEFAULT:
# DESCRIPTION: Provide a maximum obs time that shouldn't be exceeded
#
max_obs_time = ""

####### min_obs_time #######
#
# NAME: max_obs_time
# OPTIONS:
# TYPE: string
# FORMAT: YYYY-MM-DD HH:MM:SS
# DEFAULT:
# DESCRIPTION: Provide a minimum obs time that shouldn't be exceeded
#
min_obs_time = ""

####### data_url #######
#
# NAME: data URL
# OPTIONS:
# TYPE: string
# FORMAT: "/path/to/data"
# DEFAULT:
# DESCRIPTION: Provide a URL to a dataset you want to match times to
#
data_url = ""

"""
