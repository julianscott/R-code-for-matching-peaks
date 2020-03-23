# Code for:

(1) identifying peaks in two data time series, 

(2) matching peaks separated by a lag/lead for upstream/downstream timestamps, 

(3) plotting peaks.

# Potential uses of this script includes: 

(1) automatically identifying paired peaks in two dependent data time series, 

(2) calculating discharge-dependent lead or lag 'shift' between stream gage observations and upstream or downstream river stage observations, 

(3) constructing stage-discharge rating tables or rating curves at a study reach with an upstream/downstream stream gage, etc.

## SCRIPTS:

# peak_matching_script.R - 
Use this script for matching peaks in two properly formatted time series.

# functions_for_peak_matching.R - 
Download this script and place in working directory or some other noted directory. This script provides SOURCE code for functions used in peak_matching_script.R.

# sync_time_series_script.R -
Use this optional script for creating input that is properly formatted for analysis in the peak_matching_script.R.

## EXAMPLE DATA:

# example_synced_sensor_data.csv - 
Example data for use in running the peak_matching_script.R 

# example_stage_sensor_data.csv - 
Example data for use in running the sync_time_series_script.R.




                        
                            

