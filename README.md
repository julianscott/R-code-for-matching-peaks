# R-code-for-matching-peaks
# Author: Julian A Scott, National Stream and Aquatic Ecology Center, US Forest Service
# Contact: julianscott@usda.gov - 970-295-5974 - julianscotta@gmail.com
# 23 March, 2020

# Code for (1) identifying peaks in two data time series, (2) matching peaks separated by 
# a lag/lead for upstream/downstream timestamps, (3) plotting peaks.

# Potential uses of this script includes 
# (1) automatically identifying paired peaks in 
# two dependent data time series, 
# (2) calculating discharge-dependent lead or lag 'shift' 
# between stream gage observations and upstream or downstream river stage observations, 
# (3) constructing stage-discharge rating tables or rating curves at a study reach with an
# upstream/downstream stream gage, etc.

# Required input file:
# 1. A data frame (called 'pt_data' below) with three columns: datetime, stage, and q.
#    File can be a user-provided csv or produced using the 'sensor_and_gage_sync_generalized' 
#    script available at https://github.com/julianscott/R-code-for-matching-peaks.
     
#    See example data frame (provided below at the github link) to observe required format.


# To run script, either download pt_data download example synced up pressure stage and discharge 
# time series data (provided below at the github link). 
# Or, read in your own synced up data. Format of your own data must match example data. 
# pt_data <- fread("https://raw.githubusercontent.com/julianscott/R-code-for-matching-peaks/master/synced_Q_and_stage_timeseries.csv")

# This strategy uses 4 steps for matching peaks. Starting with the dataframe called pt_data,
# the dataframe name is modified after each step so that, at the end, we have objects called
# pt_data2, pt_data3, pt_data4, and pt_data5. Each consequtive dataframe is based on the 
# previous table and has the same number of rows, but additional columns, containing slope,
# peak, and other types of diagnostic data. 

# pt_data: Dataframe with two time series of data that shair a datetime column.
# pt_data2: Add moving average slope, using a defined window. Also adds a slope sign column.
# pt_data3: Flag each observation that meets the criteria as a 'peak'. Uses a rolling window
#           to tests whether the observation is the maximum in the defined window.
# pt_data4: For each flagged 'peak', flag again those observations that occur within a window where 
#           the slope is greater than the given percentile of all slopes (95% is default).  
# pt_data5: Tests whether there are flagged peaks in one time series that precede or follow the 
#           occurrence of a flagged peak in the other time series. Takes arguments to set the 
#           window in which to look for peaks and for which peak flag to base the test on: 
#           "basic peaks" from pt_data3 or "peaks that occurs within a window with a particularly 
#           steep slope", from pt_data4)


# A final object called pt_data6 is constructed that is in 'long' form, for plotting purposes, 
# with stage and discharge data stacked on top of each other, so there are twice as many rows 
# as the other pt_data* objects.  

# Finally, the objects "matched_peaks" is created. This object contains the results of this
# peak matching script, with the datetime and magnitude of the paired peaks and the shift
# in time (lead or lag) between each paired peak.

# Note the use of scaled magnitudes throughout this script. This is a key component for
# visuallizing paired peaks for scalars that can vary by orders of magnitude.
# The scale::rescales() function rescales each time series, by group or overal, to the 
# range 0 to 1. The actual magnitudes are kept throughout, though the scaled values are
# used exclusively in plotting.

# The intent is for the analyst to use the arguments in these functions to iteritively search
# for the "best" settings for a given study - ONE SETTING IS NOT BEST FOR ALL SITUATIONS.
