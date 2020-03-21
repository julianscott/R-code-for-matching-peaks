###############################################################################
# Code for identifying peaks in two timeseries of data, matching peaks separated by 
# a lag/lead for upstream/downstream timestamps, and plotting for QAQC


packages <- c("RCurl","oce","roll","RcppRoll","dataRetrieval","formattable","data.table","tibbletime","hydroTSM","xtable","tidyverse","chron","lubridate","ggpubr","segmented","stargazer","zoo","readxl")

# Check to see if each is installed, and install if not.
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))
}

# now, use the lapply function to load the installed libraries in the packages list
lapply(packages,library,character.only=TRUE)


# Input file is:
# 1. A data frame (called 'pt_data' below) with three columns: datetime, stage, and q.
#    File can be a user-provided csv or produced using the sensor_and_gage_sync_generalized script available at
#    https://github.com/julianscott/R-code-for-matching-peaks


# To run script, download pt_data download example synced up pressure stage and discharge timeseries data. 
# Or, read in your own synced up data. Format of your own data must match example data. 
pt_data <- fread("https://raw.githubusercontent.com/julianscott/R-code-for-matching-peaks/master/synced_Q_and_stage_timeseries.csv")

# Ensure datetime data is formatted correctly. For the example, its year-month-day hour:minute:second 
# in the "America/Los_Angeles" time zone. 

# Note, convieniently, lubridate can handle daylight saving time.
# ?tz
# check OlsonNames for list of valid time zones
# OlsonNames(tzdir = NULL)

# proj_tz = "America/Los_Angeles" #
proj_tz = "Etc/GMT-8" #
pt_data$DateTime <- lubridate::ymd_hms(pt_data$DateTime,tz = proj_tz)

# Define the number of chunks to split the time series into for viewing
cut_number <- 5

# Format the dataframe with some useful columns
pt_data <- pt_data %>%
  dplyr::select(DateTime,h,q_cfs) %>%
  # Add column for discharge in cubic meters/second.
  # Add a catagorical variable that splits the time series into chucks (for plotting)
  mutate(q = q_cfs/(3.2808^3),
        cut_ts = cut(DateTime,cut_number)) %>%
  # Group by cut_ts and add two columns to show h and q data rescaled to the range 0-1 (for visualizing peaks))
  group_by(cut_ts) %>%
  mutate(q_rs = scales::rescale(h),
         h_rs = scales::rescale(q))


# Graphically inspect q and h data by plotting rescaled variables over time.
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
names(cbPalette) <- as.character(expression(grey,orange,lightblue,green,yellow,darkblue,orange,purple))
model_color <- cbPalette[c("grey","orange")]
names(model_color) <-as.character(c("Stage","Discharge"))

ggplot(data = pt_data) +
  # geom_point(aes(x = DateTime, y = h_rs,color = "Stage"),size = 1) +       
  # geom_point(aes(x = DateTime, y = q_rs,color = "Discharge"),size = 1) +  
  geom_line(aes(x = DateTime, y = h_rs,color = "Stage"),size = 0.5,alpha = 1) +       
  geom_line(aes(x = DateTime, y = q_rs,color = "Discharge"),size = 0.5,alpha = 1) +  
  scale_color_manual(name = '',
                     values = model_color) +
  ylab("Rescaled q and h") +
  xlab("Date") +
  theme(panel.background = element_rect(fill = "white"),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.line = element_line(size = 0.5, linetype = "solid",
                                 colour = "black")) +
  facet_wrap(~cut_ts ,scales = "free",ncol = 1)

# Same, but adjust h_rs and q_rs to scale overall, rather than within cut_ts
pt_data %>% 
  group_by() %>% # clears grouping variable
  mutate(q_rs = scales::rescale(h),
         h_rs = scales::rescale(q)) %>% 
  ggplot() +
    # geom_point(aes(x = DateTime, y = h_rs,color = "Stage"),size = 1) +       
    # geom_point(aes(x = DateTime, y = q_rs,color = "Discharge"),size = 1) +  
    geom_line(aes(x = DateTime, y = h_rs,color = "Stage"),size = 0.5,alpha = 1) +       
    geom_line(aes(x = DateTime, y = q_rs,color = "Discharge"),size = 0.5,alpha = 1) +  
    scale_color_manual(name = '',
                       values = model_color) +
    ylab("Rescaled q and h") +
    xlab("Date") +
    theme(panel.background = element_rect(fill = "white"),
          strip.background = element_blank(),
          strip.text.x = element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.line = element_line(size = 0.5, linetype = "solid",
                                   colour = "black")) +
    facet_wrap(~cut_ts ,scales = "free",ncol = 1)


#### End section of script for reading and formating syced timeseries data 
###########################################################################################
#### Begin section for matching peaks in the two timeseries

# This strategy uses 4 steps for matching peaks. Starting with the dataframe called pt_data,
# the dataframe name is modified after each step so that, at the end, we have objects called
# pt_data2, pt_data3, pt_data4, and pt_data5. Each consequtive dataframe is based on the 
# previous table and has the same number of rows, but additional columns. 

# pt_data2 

# After the addition of a
# Use a moving average functions (hydro_timeseries_analysis.R) to get sliding means, slopes, peaks, etc from
# time series of q and h



# source("E:\\_DoD\\_Camp_Pendleton_Survey\\R_code\\Analysis\\hydro_timeseries_analysis_generalized.R")
source("hydro_timeseries_analysis_generalized.R")

# calculate a rolling average and rolling slope 

# define the window for the calculation. Assumes each data point is a 15 minute measurement.
slope_window = 3 # where 1 = 15 minutes, for 15-minute data.
pt_data2 <- pt_data %>%
  # remove grouping 
  group_by() %>% 
  as_tbl_time(index = DateTime) %>%
  mutate(h_slope = .rolling_lm(x = DateTime, y = h), # sliding elevation slope (dh/dt), for value + previous 4 measures (right aligned)
         h_msign = ifelse(h_slope > 0, "+","-" ),     # sign of the elevation slope
         q_slope = .rolling_lm(x = DateTime, y = q), # sliding q slope (dq/dt), for value + previous 4 measures (right aligned)
         q_msign = ifelse(q_slope > 0, "+","-" ))     # ssign of the q slope

# identify peaks in h and q
# For example dataset, each row is 15 minutes, so  1=15 minutes, 2=30 minutes, etc.
# 23 is 5.75 hrs
# 32 is 8 hours
# 45 is 11.25 hrs
# 95 is 23.45 hrs
# 479 os ~5 days
# 2879 is ~30 days
peak_window = 32 # i.e., maximum values in an 8 hour moving window are flagged as 'peaks'
pt_data3 <- .peak_fun(df = pt_data2,
                      measure = c("h","q"),
                      date = "DateTime",
                      peak_window = peak_window)
# head(filter(pt_data3,!is.na(DateTime)))

# Produce a plot showing hydrograph of stage and discharge, with peaks marked by the peak_fun function.
# Plotting function 'plot_data3_fun' takes pt_data3, the number of facets you want to split the time series
# into, and the peak_window as arguments. It returns both the plot and the plot data as a list.
get_plot_and_plotdata3 <- .plot_data3_fun(data = pt_data3,facet_count = 5,window = peak_window)
plot_pt_data3 <- get_plot_and_plotdata3[1]
plotdata_pt_data3 <- get_plot_and_plotdata3[2]

# print plot
plot_pt_data3

# Test if the peaks marked in 'pt_data3' occur in a window of time (test window) 
# when abs(slope) was > 95% of all slopes. The intent is to distinguish those peaks 
# that are associated with dramatic increases (i.e. true peaks)

test_window = 5
slope_percentile = 0.95
pt_data4 <- .peak_test(df = pt_data3,
                       slope = c("h_slope","q_slope"),
                       slope_percentile = slope_percentile,
                       peak = c("h_peak_logi","q_peak_logi"),
                       date = "DateTime",
                       test_window = test_window)

get_plot_and_plotdata4 <- .plot_data4_fun(data = pt_data4,
                                          facet_count = 5,
                                          window = peak_window,
                                          groupscale = T)
plot_pt_data4 <- get_plot_and_plotdata4[1]
plotdata_pt_data4 <- get_plot_and_plotdata4[2]
plot_pt_data4


# For every marked peak h, look in the match_window for a marked peak q.
# search_direction depends on whether the stream gage is upstream ("right"
# or downstream ("left") of the sensor. If a marked q is found, return the 
# number of rows separating the h peak from the q peak. This is the 'match_row' column.

# The variables q_logi and h_logi allow you to set which logic variable to base the test on:
# q_peak_logi and h_peak_logi are tests of whether the measurement is the maximum in the window 
# (i.e. a marked peak), while q_slope_test and h_slope_test are tests of whether the slopes in 
# the graph around the marked peak exceed the defined slope threshold (i.e. the latter 
# is usually more stringent the former).

# add in overrides if desired
# manual additions
# DT_pq_m1 <- ymd_hm(c("2019-02-14 15:15","2019-02-14 09:15"),tz=proj_tz)
# DT_ph_m1 <- ymd_hm(c("2019-02-14 15:30","2019-02-14 09:30"),tz=proj_tz)
# override_df <- data.frame(DT_pq = DT_pq_m1,DT_ph = DT_ph_m1)
# pt_data4 <- .override_fun(dt_df = override_df,record_df = pt_data4)

match_window = 23 
search_direction = "right" # stream gage is upstream ("right") or downstream ("left") of the sensor
pt_data5 <- .match_test(df = pt_data4,
                        # q_logi = "q_peak_logi",
                        # h_logi = "h_peak_logi",
                        q_logi = "q_slope_test",
                        h_logi = "h_slope_test",
                        date = "DateTime",
                        search_direction = search_direction, # 
                        test_window = match_window)
# when a marked peak in h does not have a matching marked peak q, -Inf is produced. Replace with NA
# Also, remove rows with DateTime NA, which are padded rows at the head and tail, required for rolling analysis
# in the function.
pt_data5 <- pt_data5 %>% 
  mutate(match_row = replace(match_row,match_row == -Inf,NA)) %>% 
  filter(!is.na(DateTime))

# Add a column that indicates, for each qualifying peak H, the row number of the qualifying
# matched peak Q.
if(search_direction == "left"){
  pt_data5 <- mutate(pt_data5,pq_row = row_number() + match_row)
} else if(search_direction == "right") {
  pt_data5 <- mutate(pt_data5, pq_row = row_number() - match_row) 
}

get_plot_and_plotdata5 <- .plot_data5_fun(data = pt_data5,
                                          facet_count = 5,
                                          window = peak_window,
                                          groupscale = T)
plot_pt_data5 <- get_plot_and_plotdata5[1]
plotdata_pt_data5 <- get_plot_and_plotdata5[2]
plot_pt_data5


matched_peaks <- .matched_peaks_fun(df = pt_data5,
                                    search_dir = search_direction,
                                    match_window = match_window,
                                    view_minutes = 23)



########################################################
# matched_peaks table contains everying needed for evaluating the marked peaks
# code below is for visualizing time series with marked peaks and plotting 
# rating curves and time-shift vs discharge curves. 

# rename for work
sensor_data <- pt_data5 %>%
  select(-c(h_slope,q_slope))

# get intervals of time for each identified peak
# create a tbl of peak number + viewing interval
peak_interval <- matched_peaks %>%
  select(peak_n,cut_peak_n,view_final)

# create a table of original DateTime (no time shift)
DT_table <- sensor_data %>%
  select(DateTime)

# cross DT_table and peak_interval to create a dataframe with length equal 
#  to nrow(DT_table) * nrow(peak_interval), with all combinations of (DateTime) 
#  and view_final. Then, test whether each DateTime is in the given view_final,
#  mark those rows where the DateTime %within% view_final == TRUE, drop those that are not.
#  This effectively identifies a viewing window for each peak.  Join the resulting table
#  to sensor_data.

sensor_data2 <- crossing(DT_table,peak_interval) %>%
  mutate(active = if_else(DateTime %within% view_final,peak_n,NA_integer_)) %>%
  drop_na(active) %>%
  select(DateTime,peak_n,cut_peak_n) %>%
  # distinct(DateTime,.keep_all = T) %>%
  right_join(sensor_data,by = "DateTime") %>%
  group_by(cut_peak_n,peak_n) %>%
  # this rescales q and h data to the PANEL
  mutate(panel_qrs= scales::rescale(q),
         panel_hrs = scales::rescale(h)) %>%
  ungroup()
# sensor_data3 <- distinct(sensor_data2,DateTime,.keep_all = T)

# add panel_qrs and panel_hrs to matched_peaks
matched_peaks <- matched_peaks %>%
  select(-view_final) %>%
  # look up panel_qrs by DateTime/DT_pq
  inner_join(rename(select(sensor_data2,DateTime,panel_qrs),DT_pq = DateTime),by = "DT_pq") %>%
  # look up panel_hrs by DateTime/DT_ph
  inner_join(rename(select(sensor_data2,DateTime,panel_hrs),DT_ph = DateTime),by = "DT_ph")


################################ make plot of all data
names(cbPalette) <- as.character(expression(grey,orange,lightblue,green,yellow,darkblue,orange,purple))
model_color <- cbPalette[c("grey","lightblue","purple","orange")]
names(model_color) <-as.character(c("Discharge","Stage","Q Peak","H Peak"))

ggplot(data = sensor_data) + 
  geom_point(aes(x = DateTime, y = q_rs,color = "Discharge"),size = 0.5) +
  geom_point(aes(x = DateTime, y = h_rs,color = "Stage"),size = 0.5) +
  geom_point(data = matched_peaks,
             aes(x = DT_pq, y = pq_rs,color = "Q Peak"),size = 1) +
  geom_point(data = matched_peaks,
             aes(x = DT_ph, y = ph_rs,color = "H Peak"),size = 1) +
  scale_color_manual(name = '',
                     values = model_color,
                     labels =c("Discharge","H Peak","Q Peak","Stage")) +
  ylab("Rescaled h and q") +
  xlab("Date") +
  # ggtitle(paste0(round(peak_window*15/60/24),"day peak /",
  #                slope_percentile,"% slope cutoff / ",
  #                round(match_window*15/60),"hr match window")) +
  theme_minimal(base_size =12) +
  # theme_void(base_size =24) +
  theme(strip.background = element_blank(), 
        strip.text.x = element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.line = element_line(size = 0.5, linetype = "solid",
                                 colour = "black")) +
  facet_wrap(~cut_ts ,scales = "free",ncol = 1)


################################ make plot of just peaks
# set plot_i and plot_f to control how many plots to print  
max(sensor_data2$peak_n,na.rm = T)
plot_i = 1
plot_f = 1
ggplot(data = filter(sensor_data2,!is.na(peak_n) & peak_n %in% plot_i:plot_f)) +
  geom_point(aes(x = DateTime, y = panel_qrs,color = "Discharge"),size = 0.5) +
  geom_point(aes(x = DateTime, y = panel_hrs,color = "Stage"),size = 0.5) +
  geom_point(data = filter(matched_peaks,peak_n %in% plot_i:plot_f),
             aes(x =  DT_pq, y = panel_qrs,color = "Q Peak"),size = 1) +
  geom_point(data = filter(matched_peaks,peak_n %in% plot_i:plot_f),
             aes(x = DT_ph, y = panel_hrs,color = "H Peak"),size = 1) +
  scale_color_manual(name = '',
                     values = model_color,
                     labels =c("Discharge","H Peak","Q Peak","Stage")) +
  ylab("Rescaled h and q") +
  # xlab("Date") +
  # ggtitle(paste0(round(peak_window*15/60),"hour peak /",
  #                slope_percentile,"% slope cutoff / ",
  #                round(match_window*15/60),"hr match window")) +
  theme_minimal(base_size =12) +
  theme(strip.background = element_blank(),
        # strip.text.x = element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.line = element_line(size = 0.5, linetype = "solid",
                                 colour = "black"))+
  facet_wrap(cut_peak_n ~ peak_n, scales = "free")

################################ make plot of rating curve
ggplot()+
  geom_point(data = sensor_data, aes(x = q,y=h,color = "raw data"))+
  geom_point(data= matched_peaks, aes(x = pq,y = ph,color = "Lag Adjusted Peaks"))+
  # ggtitle(paste0(round(peak_window*15/60),"hr peak /",
  #                slope_percentile,"% slope cutoff / ",
  #                round(match_window*15/60),"hr match window"))+
  scale_color_manual(name = '',
                     values = c("red","dark grey"),
                     labels =c("Lag Adjusted Peaks","raw data")) +
  # theme_minimal(base_size =18) +
  theme(legend.position="bottom")

################################ make plot of shift vs q plot
shift_df <- matched_peaks
# filter(!peak_n %in% c(15,19,20))

ggplot()+
  geom_point(data = shift_df, aes(x = pq,y = shift),color = "red")+
  theme(legend.position="bottom")


################################ filter to inspect time periods
sens_filter2 <- filter(sensor_data2,between(DateTime,ymd("2018-11-29",tz=proj_tz),ymd("2018-12-01",tz=proj_tz)))
match_filter <- filter(matched_peaks,between(DT_pq,ymd("2018-11-29",tz=proj_tz),ymd("2018-12-01",tz=proj_tz)))

# sens_filter2_date <- filter(sensor_data2,peak_n == 12)
sens_filter2 <- filter(sensor_data2,peak_n == 10)
# View(sens_filter2)
match_filter <- filter(matched_peaks2,peak_n == 10)
# View(match_filter)

ggplot(data = sens_filter2) + 
  geom_point(aes(x = DateTime, y = q_rs),color = "black",size = 0.5) +
  geom_point(aes(x = DateTime, y = h_rs),color = "light grey",size = 0.5) +
  geom_point(data = match_filter,
             aes(x = DT_pq, y = pq_rs),color = "green",size = 1) +
  geom_point(data = match_filter,
             aes(x = DT_ph, y = ph_rs),color = "red",size = 1) +
  ylab("Rescaled h and q") +
  # xlab("Date") +
  # ggtitle(paste0(round(peak_window*15/60/24),"day peak /",
  #                slope_percentile,"% slope cutoff / ",
  #                round(match_window*15/60),"hr match window")) +
  theme_minimal(base_size =12) +
  # theme_void(base_size =24) +
  theme(strip.background = element_blank(), 
        strip.text.x = element_blank(),
        # panel.grid.major=element_blank(),
        # panel.grid.minor=element_blank(),
        axis.line = element_line(size = 0.5, linetype = "solid",
                                 colour = "black"))


