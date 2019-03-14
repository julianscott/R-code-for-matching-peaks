packages <- c("oce","roll","RcppRoll","dataRetrieval","formattable","data.table","tibbletime","hydroTSM","xtable","tidyverse","chron","lubridate","ggpubr","segmented","stargazer","zoo","readxl")

# Check to see if each is installed, and install if not.
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))
}

# now, use the lapply function to load the installed libraries in the packages list
lapply(packages,library,character.only=TRUE)

####################################################
#### read in PT data 
#### PT data should have two columns, datetime and stage, in that order
####################################################
# choose.files()
# setwd("E:\\_DoD\\_Camp_Pendleton_Survey\\Hydrology\\Pressure_Transducer_data\\PT_data_Time_Flow_Alignment_work\\031219_fullrecord\\")
setwd("E:\\_DoD\\_Camp_Pendleton_Survey\\R_code\\Analysis\\Beasley")
PTfiles <- list.files(pattern = "*.csv")
PTfiles
# length(PTfiles)

# PT stream sensor files should have two columns, 1=datetime,2=stage
PT <- read.csv(PTfiles[[2]],header=T,stringsAsFactors = FALSE)
colnames(PT) 

# standardize column names
colnames(PT) <- c("DateTime","h")

# proj_tz = "America/Los_Angeles" # 
proj_tz = "America/Phoenix" # 

# Format datetimes. Check input date data - is it mdy_hm or ymd_hm? change to match
head(PT)
PT <- mutate(PT,DateTime = mdy_hm(DateTime,tz = proj_tz))

# view interval of time series
# interval(min(ymd_hms(PT[,"DateTime"],tz = proj_tz)),max(ymd_hms(PT[,"DateTime"],tz = proj_tz)))
PT_start <- min(ymd_hms(PT[,"DateTime"],tz = proj_tz))
PT_end <- max(ymd_hms(PT[,"DateTime"],tz = proj_tz))

# filter PT data to a certain date range if desired
# PT_start <- mdy_hm("2018-03-11 11:00",tz = proj_tz)
# record max to use 2018-12-03 22:45:00 datetime8
# PT_end <- ymd_hm("2018-12-03 22:45",tz = proj_tz)
# PT <- PT %>%
#   filter(between(DateTime,PT_start,PT_end))

unique(minute(PT$DateTime))
unique(second(PT$DateTime))
unique(hour(PT$DateTime))

# fill in NA values if appropriate/needed
PT[which(is.na(PT$h)),]  
PT[which(is.na(PT$DateTime)),]  

# fix erroneous measurements if needed
# filter(PT_2,DateTime ==  ymd_hm("2018-04-22 18:45",tz = proj_tz))$h_2
# PT_2[PT_2$DateTime ==  ymd_hm("2018-04-22 18:00",tz = proj_tz),]$h_2 <- 77.132

# rename for work
pt_data <- PT

#### Code for reading in 15 minuted data from usgs gages

# setwd("E:\\_DoD\\_Camp_Pendleton_Survey\\Hydrology\\WSE_data\\asof032018")

# read in flow data direct from USGS website
# readNWISuv acquires the current/historical observations (15 minute data)

# tz7 = "America/Denver" # GTM-7 (datetime7) ?tz
# proj_tz = "America/Denver" # 
# proj_tz = "America/Phoenix" # 
# ?tz
# check OlsonNames for list of valid time zones
# OlsonNames(tzdir = NULL)
# proj_tz

# stream_gage_start_date = "2017-06-28 0:00"
# stream_gage_end_date = "2018-10-11 23:45"
# interval(min(ymd_hms(PT[,"DateTime"],tz = proj_tz)),max(ymd_hms(PT[,"DateTime"],tz = proj_tz)))

# 11044300 is fallbrook
# 09506000 is beasley
gage1_raw <- readNWISuv(site = '09506000', parameterCd="00060",
                        startDate = date(PT_start),
                        endDate = date(PT_end),
                        tz = proj_tz)
 
gage2_raw <- readNWISuv(site='11044350', parameterCd="00060",
                        startDate = date(PT_start),
                        endDate = date(PT_end),
                        tz = proj_tz)

# I am adding two gages together at the SMR, so I need to align the timestamps first
# If you don't need to do this, just do A.

# Rename columns and organize for analysis 
# Be certain of your time zones!! This script assumes all data is in the same time zone, proj_tz

# A.
gage1 <- dplyr::select(gage1_raw,site_no,dateTime,X_00060_00000,tz_cd)
colnames(gage1) <- c("site_n","DateTime","gage1_cfs","tz") 

# B.
gage2 <- dplyr::select(gage2_raw,site_no,dateTime,X_00060_00000,tz_cd)
colnames(gage2) <- c("site_n","DateTime","gage2_cfs","tz") 

# to ensure a continuous sequence of dates, from start to end, 
# by 15 minute interval, I create my own sequence. Importantly, 
# this method accounts for daylight savings time (e.g. observe behaviour
# in mid march and early november)
dateseq <- seq(min(gage1$DateTime),
               max(gage1$DateTime), 
               by = '15 mins')

# create new dataframe, with my dateseq as the first column
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
# Use this code for processing 1 usgs gage
###########################

hydf <- hydf_dates %>%
  dplyr::left_join(select(gage1,DateTime,gage1_cfs),by = "DateTime") %>%
  mutate(q_cfs = gage1_cfs) %>%
  select(DateTime,q_cfs)
head(hydf)

#### End section of script for reading in a formating gage data and pressure transducer data
###########################################################################################
#### Begin PT linkage to gage data

q_date <- hydf$DateTime  # date record of Q
q <- hydf$q_cfs      # gage record of Q
# use linear interpolation (the approx() command) to calulate the gage Q for each h in pt_data
# pt_q is a vector the same length as pt_data$datetime, with 
#  one interpolated Q for every datetime in pt_data$datetime
#  rule = 1 provides an NA for any pt_data$DateTime that is out of the q_date range
pt_q <- approx(q_date,q,xout = pt_data$DateTime,rule = 1)$y
# approx(q_date,q,xout = pt_data$DateTime[2],ties = mean,yright = 9999,yleft = -9999)
pt_data$q_cfs <- pt_q

# This is the end of the R method for syncing datetime stamps for PT and gage data
###############################################################################