
# Daylight savings tool
# Hydrologic sensors, like pressure transducers by Onset, usually do not switch from standard to daylight savings time. So, 
# for example, when a sensor is deployed on March 10th, 2018, and runs until April, 2018, the daylight savings time 
# change that occurs at 2 am on March 11th is not noted. Rather, the sensor maintains whichever timezone hour offset
# from GMT specified (GMT-8, for Pacific Standard Time, for example) for all measurements. Best practice is to download stream 
# discharge data from a USGS gage and specify the same time zone as the sensor, however, if there is need to change from
# a standard (non-daylight savings) time zone, into a daylight savings time zone (or visa versa), here is how you can do it.

packages <- c("oce","roll","RcppRoll","dataRetrieval","formattable","data.table","tibbletime","hydroTSM","xtable","tidyverse","chron","lubridate","ggpubr","segmented","stargazer","zoo","readxl")

# Check to see if each is installed, and install if not.
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))
}

# now, use the lapply function to load the installed libraries in the packages list
lapply(packages,library,character.only=TRUE)



# Sensor data is 8 hours behind GMT (confusingly, this tz is 'Etc/GMT+8', the sign is intentionally inverted)
# see: https://en.wikipedia.org/wiki/Tz_database#Area
pt_data <- fread("https://raw.githubusercontent.com/julianscott/R-code-for-matching-peaks/master/synced_Q_and_stage_timeseries.csv")
original_tz = "Etc/GMT+8" # Does not switch over for daylight savings time
new_tz = "America/Los_Angeles" # Does swtich over for daylight savings time

# Format DateTime to known original timezone
pt_data$DateTime <- lubridate::ymd_hms(pt_data$DateTime,tz = original_tz)

# Convert to timezone that considers daylight savings
pt_data$DateTime_new <- format(pt_data$DateTime, tz=new_tz,usetz=TRUE)

# Convert from character to date
pt_data$DateTime_new <- lubridate::ymd_hms(pt_data$DateTime_new,tz = new_tz)
str(pt_data)

# Convert back to original time zone
# pt_data$DateTime_original <- with_tz(pt_data$DateTime_new, tzone=original_tz)
# pt_data$DateTime_original <- lubridate::ymd_hms(pt_data$DateTime_new,tz = original_tz)
# pt_data$DateTime_original <- format(pt_data$DateTime, tz=original_tz,usetz=TRUE)
# pt_data$DateTime_original <- as.POSIXct(format(pt_data$DateTime_new, tz=original_tz), tz=original_tz)
pt_data$DateTime_original <- strptime(pt_data$DateTime_new,format = )



# Observe daylight savings behavior
hydf_dates %>% 
# pt_data %>%
  filter(between(DateTime,
                 ymd_hm("2018-11-04 00:00",tz = original_tz),
                 ymd_hm("2018-11-04 04:00",tz = original_tz)))
str(pt_data$DateTime)
summary(pt_data$DateTime)
strptime()
tmp

tmp$DateTime2 <-  format(tmp$DateTime, tz="America/Los_Angeles",usetz=TRUE)
tmp$DateTime3 <-  format(tmp$DateTime, tz="Etc/GMT+8",usetz=TRUE)
tmp$DateTime4 <-  force_tz(tmp$DateTime, tz="America/Los_Angeles")
# format(pb.date,usetz=TRUE, tz="Etc/GMT+8")
tmp
str(tmp)

nowAsCharacter = format(now, tz="GMT")
as.POSIXct(format(nowAsCharacter, tz="GMT"), tz="GMT")


now = Sys.time()
nowAsCharacter = format(now, tz="GMT")
as.POSIXct(format(nowAsCharacter, tz="GMT"), tz="GMT")



df <- data.frame(actualtime = c(
  '2015-04-15 13:10:00',
  '2015-04-15 14:22:00',
  '2015-04-15 10:14:00'),
  timezone = c(
    'Australia/Sydney',
    'Australia/Perth',
    'Australia/Perth'))

raw <- as.POSIXct(strptime(
  df$actualtime[2],
  format = "%Y-%m-%d %H:%M:%S", 
  tz ="Australia/Sydney"),
  tz = "Australia/Sydney")

df
raw
# starts conversion
converted <- as.POSIXlt(raw, tz = df$timezone[2])


# From Tobias: https://stackoverflow.com/questions/13865172/handling-dates-when-we-switch-to-daylight-savings-time-and-back
# Convert datetimes to account for daylight savings 
# as.POSIXct.no.dst <- function (x, tz = "", format="%m/%d/%Y %H:%M", offset="+0100", ...)
# {
#   x <- paste(x, offset)
#   format <- paste(format, "%z")
#   as.POSIXct(x, tz, format=format)
# }
# 
# DateTime_MDT <- as.POSIXct.no.dst(x = ld_q_raw$DateTime,tz = local_tz,offset = "-0700")
# ld_q_raw$DateTime_MDT <- DateTime_MDT



