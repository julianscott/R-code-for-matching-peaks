###############################################################################
# This code can be used to sync up time series of data. Specifically, it is for
# syncing time series of pressure transducer stage data  with upstream or 
# downstream stream flow gage data. An option is included to combine flow from 
# a mainstem gage with a nearby tributaries into a single flow record, prior
# to syncing with the sensor data.

# The formatted results of this code can be used directly in the peak_matching_script.R
# and functions_for_peak_matching.R files found at: 
# https://github.com/julianscott/R-code-for-matching-peaks


# Script objectives
# 1. read and format pressure transducer data 
# 2. download and format usgs stream gage data from 1 or 2 gages.
# 3. sync the date and time for all time series of data

packages <- c("dataRetrieval","data.table","tidyverse","lubridate")

# Check to see if each is installed, and install if not.
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))
}

# now, use the lapply function to load the installed libraries in the packages list
lapply(packages,library,character.only=TRUE)

####################################################
#### read in pressure transducer (PT) data 
#### PT data should have two columns, datetime and stage, in that order
####################################################

# Working directory taken from directory of the R-code-for-matching-peaks-master.Rproj file.
# Else, use setwd.

# read in example pressure transducer data from my github site:
# Like this example, PT sensor files should have two columns, 1=datetime,2=stage
PT <- fread("https://raw.githubusercontent.com/julianscott/R-code-for-matching-peaks/master/example_stage_sensor_data.csv")

# standardize column names
colnames(PT) <- c("DateTime","h")

# Mangement of the date and time data is absolutely critical for time series analysis.
# Daylight savings must be considered (the R package 'lubridate' is great for this).

# Set time zone for the project. 
proj_tz = "Etc/GMT+8" # This does not adjust for daylight savings, which matches how most hydrologic
# sensors collect data.

# ?tz
# check OlsonNames for list of valid time zones
# OlsonNames(tzdir = NULL)

# Format datetimes. Check input date data - is it mdy_hm or ymd_hm? change to match
PT <- mutate(PT,DateTime = mdy_hm(DateTime,tz = proj_tz))

# view interval of time series
interval(min(ymd_hms(PT[,"DateTime"],tz = proj_tz)),max(ymd_hms(PT[,"DateTime"],tz = proj_tz)))

# We will limit our analysis to the interval defined by the start and stop of the PT data
PT_start <- min(ymd_hms(PT[,"DateTime"],tz = proj_tz))
PT_end <- max(ymd_hms(PT[,"DateTime"],tz = proj_tz))

# QAQC tip - view the unique second, minute, and hours that are in your dataset. Is this as you expected?
# Datetimes and formats are finicky. This dataset uses 15 minute intervals.
unique(second(PT$DateTime))
unique(minute(PT$DateTime))
unique(hour(PT$DateTime))

# QAQC tip - are there gaps in the sensor time series? Here, I find 45 minutes of stage data missing. 
# I fill in the data with the mean of the bounding data (77.047). 

#  check for rows with no stage measurements
PT[which(is.na(PT$h)),]  
# check for rows with no datetime 
PT[which(is.na(PT$DateTime)),]  

# fix erroneous measurements if needed
PT[8356:8360,]
PT[PT$DateTime ==  ymd_hm("2018-06-06 08:30",tz = proj_tz),"h"] <- 77.047
PT[PT$DateTime ==  ymd_hm("2018-06-06 08:45",tz = proj_tz),"h"] <- 77.047
PT[PT$DateTime ==  ymd_hm("2018-06-06 09:00",tz = proj_tz),"h"] <- 77.047

# rename dataset for work
pt_data <- PT

#### Code for reading in 15 minuted data from usgs gages

# read in flow data for two gages direct from USGS website
# readNWISuv acquires the current/historical observations (15 minute data)
gage1_raw <- readNWISuv(site = '11044300', parameterCd="00060",
                        startDate = date(PT_start),
                        endDate = date(PT_end),
                        tz = proj_tz)

gage2_raw <- readNWISuv(site='11044350', parameterCd="00060",
                        startDate = date(PT_start),
                        endDate = date(PT_end),
                        tz = proj_tz)


# Rename columns and organize for analysis.  Be certain of your time zones!! 
# This script assumes all data is in the same time zone, proj_tz.
# For this example, we are adding the flow from the mainstem (gage1) to a nearby tributary (gage2),
# so I need to align the timestamps first. If you don't need to do this, just do A.

# A.
gage1 <- dplyr::select(gage1_raw,site_no,dateTime,X_00060_00000,tz_cd)
colnames(gage1) <- c("site_n","DateTime","gage1_cfs","tz") 

# B.
gage2 <- dplyr::select(gage2_raw,site_no,dateTime,X_00060_00000,tz_cd)
colnames(gage2) <- c("site_n","DateTime","gage2_cfs","tz") 


# To ensure a continuous sequence of dates, from start to end, by 15 minute interval, I create my own sequence. 
# Importantly, this method accounts for daylight savings time. Observe behaviour on Sunday March 11th 2018 at 
# 2 am and Sunday, Nov 4th 2018 at 2 am.
dateseq <- seq(min(gage1$DateTime),
               max(gage1$DateTime), 
               by = '15 mins')

# create new dataframe, with dateseq as the first column
hydf_dates <- data.frame(DateTime = dateseq)

###########################
# Use this code for processing 2 usgs gages
###########################

# Left_join to lookup each value in the two gage records by datetime8.
hydf <- hydf_dates %>%
  dplyr::left_join(select(gage1,DateTime,gage1_cfs),by = "DateTime") %>%
  dplyr::left_join(select(gage2,DateTime,gage2_cfs),by = "DateTime") %>%
  mutate(q_cfs = gage1_cfs + gage2_cfs,
         tz_cd = dst(DateTime)) %>%
  select(DateTime,gage1_cfs,gage2_cfs,q_cfs)
head(hydf)

###########################
# Uncomment and use this code for processing 1 usgs gage
###########################

# hydf <- hydf_dates %>%
#   dplyr::left_join(select(gage1,DateTime,gage1_cfs),by = "DateTime") %>%
#   mutate(q_cfs = gage1_cfs) %>%
#   select(DateTime,q_cfs)
# head(hydf)

#### End section of script for reading and formating gage data and pressure transducer data
###########################################################################################
#### Begin syncing PT sensor and gage timeseries data

# vector of dates from Q time series
q_date <- hydf$DateTime 

# vector of flow from Q time series
q <- hydf$q_cfs      

# use linear interpolation (the approx() command) to calculate the gage Q for each h in pt_data
# pt_q is a vector the same length as pt_data$datetime, with one interpolated Q for every 
# datetime in pt_data$datetime.

# rule = 1 provides an NA for any pt_data$DateTime that is out of the q_date range
pt_q <- approx(q_date,q,xout = pt_data$DateTime,rule = 1)$y

# add vector of interpolated discahrges to the sensor time series
pt_data$q_cfs <- pt_q
# write.csv(pt_data,"synced_Q_and_stage_timeseries.csv")

# This is the end of the R method for syncing datetime stamps for PT and gage data
###############################################################################
