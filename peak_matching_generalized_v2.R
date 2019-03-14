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

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
names(cbPalette) <- as.character(expression(grey,orange,lightblue,green,yellow,darkblue,orange,purple))
model_color <- cbPalette[c("darkblue","lightblue")]
names(model_color) <-as.character(c("Stage","Discharge"))

ggplot(data = tmp) +
  geom_point(aes(x = DateTime, y = scales::rescale(h),color = "Stage")) +       # rescale standardizes both datasets to 0:1.
  geom_point(aes(x = DateTime, y = scales::rescale(q)*1.5,color = "Discharge")) +  # amplify by multiplication
  scale_color_manual(name = '',
                     values = model_color,
                     labels = c('Stage','Discharge')) +
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
search_direction = "left"
pt_data5 <- .match_test(df = pt_data4,
                        # q_logi = "q_peak_logi",
                        # h_logi = "h_peak_logi",
                        q_logi = "q_slope_test",
                        h_logi = "h_slope_test",
                        date = "DateTime",
                        search_direction = search_direction, # 
                        test_window = match_window)

source("E:\\_DoD\\_Camp_Pendleton_Survey\\R_code\\Analysis\\hydro_timeseries_analysis_generalized_complete.R")

pt_data5 <- .match_test(df = pt_data4,search_direction = search_direction)
matched_peaks <- pt_data5[[1]]
sensor_dat <- pt_data5[[2]]
sensor_dat2 <- pt_data5[[3]]

################################ make plot of all data
names(cbPalette) <- as.character(expression(grey,orange,lightblue,green,yellow,darkblue,orange,purple))
model_color <- cbPalette[c("grey","lightblue","purple","orange")]
names(model_color) <-as.character(c("Discharge","Stage","Q Peak","H Peak"))

ggplot(data = sensor_dat) + 
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
plot_i = 1
plot_f = 10
ggplot(data = filter(sensor_data2,!is.na(peak_n) & peak_n %in% plot_i:plot_f)) +
  geom_point(aes(x = DateTime, y = panel_qrs,color = "Discharge"),size = 0.5) +
  geom_point(aes(x = DateTime, y = panel_hrs,color = "Stage"),size = 0.5) +
  geom_point(data = filter(matched_peaks2,peak_n %in% plot_i:plot_f),
             aes(x =  DT_pq, y = panel_qrs,color = "Q Peak"),size = 1) +
  geom_point(data = filter(matched_peaks2,peak_n %in% plot_i:plot_f),
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
shift_df <- matched_peaks %>%
  filter(!peak_n %in% c(15,19,20))

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



