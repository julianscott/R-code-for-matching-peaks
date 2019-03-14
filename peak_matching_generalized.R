###############################################################################
# Start code for identifying peaks, adjusting peaks for lag/lead for upstream/downstream timestamps, 
# and plotting for QAQC


# To start this script, the input files are
# 1. pt_data from the results of sensor_and_gage_sync_1sensor_general.R
#  OR csv with three columns, datetime, stage, and q

# Add a auxillary columns like log value, change units, etc
pt_data <- pt_data %>%
  dplyr::select(DateTime,h,q_cfs) %>%
  mutate(q = q_cfs/(3.2808^3))

# inspect q and h data before proceding
tmp <- pt_data %>%
  mutate(cut_ts = cut(DateTime,4))

ggplot(data = tmp) +
  geom_point(aes(x = DateTime, y = scales::rescale(h)),color = "grey") +       # rescale standardizes both datasets to 0:1.
  geom_point(aes(x = DateTime, y = scales::rescale(q)*1.5),color = "black") +  # amplify by multiplication
  ylab("Rescaled q and h") +
  xlab("Date") +
  # theme_void(base_size =24)
  facet_wrap(~cut_ts ,scales = "free",ncol = 1)


# Use a moving average functions (hydro_timeseries_analysis.R) to get sliding means, slopes, peaks, etc from
# time series of q and h

#23 is 5.75 hrs
# 32 is 8 hours
# 45 is 11.25 hrs
# 95 is 23.45 hrs
# 479 os ~5 days
# 2879 is ~30 days

source("E:\\_DoD\\_Camp_Pendleton_Survey\\R_code\\Analysis\\hydro_timeseries_analysis_generalized.R")

# calculate a rolling average and rolling slope 
# define the window for the calculation
slope_window = 3 # where 1 = 15 minutes
pt_data2 <- pt_data %>%
  as_tbl_time(index = DateTime) %>%
  mutate(#h_mean = .sliding_mean(h,slope_window),      # sliding elevation mean, for value + previous 4 measures (right aligned)
         h_slope = .rolling_lm(x = DateTime, y = h), # sliding elevation slope (dh/dt), for value + previous 4 measures (right aligned)
         h_msign = ifelse(h_slope > 0, "+","-" ),     # sign of the elevation slope
         # q_mean = .sliding_mean(q,slope_window),      # sliding q mean , for value + previous 4 measures (right aligned)
         q_slope = .rolling_lm(x = DateTime, y = q), # sliding q slope (dq/dt), for value + previous 4 measures (right aligned)
         q_msign = ifelse(q_slope > 0, "+","-" ))     # ssign of the q slope

# identify peaks in h and q
peak_window = 32
pt_data3 <- .peak_fun(df = pt_data2,
                      measure = c("h","q"),
                      date = "DateTime",
                      peak_window = peak_window)
head(filter(pt_data3,!is.na(DateTime)))

# tmp <- filter(pt_data4,between(DateTime,ymd("2017-08-05",tz=proj_tz),ymd("2017-08-07",tz=proj_tz)))
# View(select(tmp,DateTime,h_peak_logi,q_peak_logi))
# View(select(tmp,DateTime,h_peak_logi,q_peak_logi,h_slope_test,q_slope_test,match_row23))
# View(tmp)
# test if peaks occur in a window of time when abs(slope) was > 95% of all slopes (i.e. true peaks)
test_window = 5
slope_percentile = 0.95
pt_data4 <- .peak_test(df = pt_data3,
                       slope = c("h_slope","q_slope"),
                       slope_percentile = slope_percentile,
                       peak = c("h_peak_logi","q_peak_logi"),
                       date = "DateTime",
                       test_window = test_window)
# for every marked peak h, look in the match_window for a marked peak q.
# search_direction depends on whether the stream gage is upstream ("right"
# or downstream ("left") of the sensor. If a marked q is found, return the 
# number of rows separating the h peak from the q peak. 

# q_logi and h_logi allow you to set witch logic variable to base the test on:
# q_peak_logi and h_peak_logi are tests of whether the measurement is the maximum in the window 
# (i.e. a marked peak), while q_slope_test and h_slope_test are tests of whether the slopes in 
# the hydro/stagegraph around the marked peak exceed the defined slope threshold (i.e. the latter 
# is more stringent the former).

# add in overrides if desired
# manual additions
DT_pq_m1 <- ymd_hm(c("2019-02-14 15:15","2019-02-14 09:15"),tz=proj_tz)
DT_ph_m1 <- ymd_hm(c("2019-02-14 15:30","2019-02-14 09:30"),tz=proj_tz)
override_df <- data.frame(DT_pq = DT_pq_m1,DT_ph = DT_ph_m1)
pt_data4 <- .override_fun(dt_df = override_df,record_df = pt_data4)

match_window = 23 
search_direction = "right"
pt_data5 <- .match_test(df = pt_data4,
                        # q_logi = "q_peak_logi",
                        # h_logi = "h_peak_logi",
                        q_logi = "q_slope_test",
                        h_logi = "h_slope_test",
                        date = "DateTime",
                        search_direction = search_direction, # 
                        test_window = match_window)

# when a marked peak in h does not have a matching marked peak q, -Inf is produced. Replace with NA
pt_data5 <- mutate(pt_data5,match_row23 = replace(match_row23,match_row23 == -Inf,NA))

# View(filter(pt_data5,date(DateTime) == "2019-02-14"))

# remove heading and tailing padded NA values used during rolling average calcs
pt_data_F <- pt_data5[!is.na(pt_data5[,"DateTime"]),]

# Depending on whether the search direction is left (sensor lags the downstream stream gage) 
# or right (sensor leads the upstream stream gage)), this code produces three new columns:
# 1. rnum is the row number of the datetime vector that will stay fixed to the stage record
# 2. sdir is a confirmation indicator of the search direction that is chosen
# 3. pq_row is the row number of a matched peak q (found within the match window, given the search direction)

if(search_direction == "left"){
  pt_data_F <- pt_data_F %>%
    mutate(rnum = row_number(),
           sdir = search_direction,
           pq_row = if_else(!is.na(match_row23),
                            rnum + match_row23 -1L,
                            NA_real_))
} else {
  pt_data_F <- pt_data_F %>%
    mutate(rnum = row_number(),
           sdir = search_direction,
           pq_row = if_else(!is.na(match_row23),
                            rnum - as.integer(match_window-match_row23),
                            NA_integer_)) 
}

# scale q and h for comparing peaks and cut time series for viewing
pt_data_F <- pt_data_F %>%
  mutate_at(vars(c(h,q)),funs(rs = scales::rescale)) %>%
  mutate(cut_ts = cut(DateTime,6))

# replace slope sign with peak designation in pt_data_F 
pt_data_F$h_msign[pt_data_F$h_peak_logi == TRUE] <- paste0(round(peak_window*15/60/24)," Day Peak")
pt_data_F$q_msign[pt_data_F$q_peak_logi == TRUE] <-  paste0(round(peak_window*15/60/24)," Day Peak")

# get rows with peak discharge measurments
pq_rows <- pt_data_F[which(!is.na(pt_data_F$pq_row)),]$pq_row

# build dataframe of just q peaks; 
# rename the following columns to indicate marked peak status:
# DateTime, q, q_rs, and rnum. Note that the rnum column is changed to
# pq_row and is the basis for the inner_join with the pt_data
q_peaks <- select(pt_data_F,rnum,DateTime,q_slope_test,q,q_rs) %>%
  filter(rnum %in% pq_rows) %>%
  rename(pq_row = rnum,DT_pq = DateTime,pq = q,pq_rs = q_rs) %>%
  select(pq_row,DT_pq,pq,pq_rs)
  
# build dataframe of h peaks and q peaks, with two date columns, one for the h and one for the lag/lead q.
# this is accomplished using the inner_join with the q_peaks table based on the pq_row column.
# Rename stage columns to indicate peak status.
matched_peaks <- filter(pt_data_F,!is.na(pq_row)) %>%
  rename(DT_ph = DateTime,ph = h,ph_rs = h_rs) %>%
  inner_join(q_peaks,by = "pq_row") %>%
  # add peak number, which can be used for idetifying a particular peak
  # add time shift between the two peaks
  mutate(peak_n = row_number(),
         cut_peak_n = cut(peak_n,2),
         shift =  as.numeric(DT_ph - DT_pq,units = "hours")) %>%
  select(peak_n,DT_pq,DT_ph,shift,pq,ph,pq_rs,ph_rs,cut_ts,cut_peak_n) %>%
  # create columns defining a period of hours before and after the marked peak to enhance QAQC visualization
  # rowwise() %>%
  mutate(target_start = pmin(DT_ph,DT_pq),
         target_end = pmax(DT_ph,DT_pq),
         view_start = target_start - match_window*minutes(15), # should not be more than match window
         view_end = target_end + match_window*minutes(15),
         view_startlag = lag(view_start,1),
         view_endlag = lag(view_end,1),
         view_startlead = lead(view_start,1),
         view_endlead = lead(view_end,1),
         view_interval = interval(start = view_start, end = view_end),
         view_interval_lag = interval(start = view_startlag, end = view_endlag),
         view_interval_lead = interval(start = view_startlead, end = view_endlead),
         int_lag = int_overlaps(view_interval,view_interval_lag),
         int_lead = int_overlaps(view_interval,view_interval_lead),
         view_start2 = if_else(int_lag == T,view_startlag,view_start),
         view_end2 = if_else(int_lead == T,view_endlead,view_end),
         view_final = interval(view_start2,view_end2),
         view_final = if_else(is.na(view_final),view_interval,view_final)) %>%
  select(-c(target_start:view_end2))
# View(matched_peaks)

# str(matched_peaks)
# all(matched_peaks$leadtmp==matched_peaks$view_lead)
# View(select(matched_peaks,peak_n,DT_ph,DT_pq,shift,view_lag,view_lead,view_final))


########################################################
# matched_peaks table contains everying needed for evaluating the marked peaks
# code below is for visualizing time series with marked peaks and plotting 
# rating curves and time-shift vs discharge curves. 

# rename for work
sensor_data <- pt_data_F %>%
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
#  This effectively identifies which view_final each peak is within.  Join the resulting table
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
matched_peaks2 <- matched_peaks %>%
  select(-view_final) %>%
  # look up panel_qrs by DateTime/DT_pq
  inner_join(rename(select(sensor_data2,DateTime,panel_qrs),DT_pq = DateTime),by = "DT_pq") %>%
  # look up panel_hrs by DateTime/DT_ph
  inner_join(rename(select(sensor_data2,DateTime,panel_hrs),DT_ph = DateTime),by = "DT_ph")


# sens_filter <- filter(sensor_data,between(DateTime,ymd("2018-03-10",tz=proj_tz),ymd("2018-03-20",tz=proj_tz)))
# match_filter <- filter(select(matched_peaks,-view_final),between(DT_pq,ymd("2018-03-10",tz=proj_tz),ymd("2018-03-20",tz=proj_tz)))

################################ make plot of all data
ggplot(data = sensor_data) + 
  geom_point(aes(x = DateTime, y = q_rs),color = "black",size = 0.5) +
  geom_point(aes(x = DateTime, y = h_rs),color = "light grey",size = 0.5) +
  geom_point(data = matched_peaks,
             aes(x = DT_pq, y = pq_rs),color = "green",size = 1) +
  geom_point(data = matched_peaks,
             aes(x = DT_ph, y = ph_rs),color = "red",size = 1) +
  ylab("Rescaled h and q") +
  # xlab("Date") +
  ggtitle(paste0(round(peak_window*15/60/24),"day peak /",
                 slope_percentile,"% slope cutoff / ",
                 round(match_window*15/60),"hr match window")) +
  theme_minimal(base_size =12) +
  # theme_void(base_size =24) +
  theme(strip.background = element_blank(), 
        strip.text.x = element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.line = element_line(size = 0.5, linetype = "solid",
                                 colour = "black"))+
  facet_wrap(~cut_ts ,scales = "free",ncol = 1)


################################ make plot of just peaks
ggplot(data = filter(sensor_data2,!is.na(peak_n))) +
  geom_point(aes(x = DateTime, y = panel_qrs),color = "black",size = 0.5) +
  geom_point(aes(x = DateTime, y = panel_hrs),color = "light grey",size = 0.5) +
  geom_point(data = matched_peaks2,
             aes(x =  DT_pq, y = panel_qrs),color = "green",size = 1) +
  geom_point(data = matched_peaks2,
             aes(x = DT_ph, y = panel_hrs),color = "red",size = 1) +
  ylab("Rescaled h and q") +
  # xlab("Date") +
  ggtitle(paste0(round(peak_window*15/60),"hour peak /",
                 slope_percentile,"% slope cutoff / ",
                 round(match_window*15/60),"hr match window")) +
  theme_minimal(base_size =12) +
  theme(strip.background = element_blank(),
        # strip.text.x = element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.line = element_line(size = 0.5, linetype = "solid",
                                 colour = "black"))+
  facet_wrap(cut_peak_n~peak_n ,scales = "free")

################################ make plot of rating curve
ggplot()+
  geom_point(data = sensor_data, aes(x = q,y=h),color = "dark grey")+
  geom_point(data= matched_peaks, aes(x = pq,y = ph),color = "red")+
  ggtitle(paste0(round(peak_window*15/60),"hr peak /",
                 slope_percentile,"% slope cutoff / ",
                 round(match_window*15/60),"hr match window"))+
  scale_color_manual(name = "", values = c("Lag Adjusted Algorithm Selected Peaks" = "red","raw data" = "dark grey"))+
  # theme_minimal(base_size =18) +
  theme(legend.position="bottom")

################################ make plot of shift vs q plot
shift_df <- matched_peaks %>%
  select(-view_final) %>%
  filter(!peak_n %in% c(15,19,20))
ggplot()+
  geom_point(data = shift_df, aes(x = pq,y = shift),color = "red")+
  theme(legend.position="bottom")


################################ filter to inspect time periods
sens_filter2 <- filter(sensor_data2,between(DateTime,ymd("2018-11-29",tz=proj_tz),ymd("2018-12-01",tz=proj_tz)))
match_filter <- filter(select(matched_peaks,-view_final),between(DT_pq,ymd("2018-11-29",tz=proj_tz),ymd("2018-12-01",tz=proj_tz)))

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



