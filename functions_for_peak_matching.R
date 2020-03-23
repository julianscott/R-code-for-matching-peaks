# hydrologic record time series analysis


.sliding_mean <- function(measure,window){
  y = rollapply(measure,width=window,FUN = mean,fill = NA,align = "right")
  return(y)
}


.rolling_lm <- rollify(.f = function(x = DateTime,y) 
  {
  tryCatch({coef(lm(y ~ x))[2]},error = function(e){NA})
}, 
window = slope_window, na_value = NA,unlist = TRUE)


# df <- pt_data2
# measure = c("h","q")
# peak_window = 479
# date = "DateTime"

.pt_data3_fun <- function(df,measure = c("h","q"),peak_window = 32,date = "DateTime"){
  measure1 = measure[1]
  head_pad <- df[0,]
  head_pad[nrow(head_pad)+round(peak_window/2),] <- NA
  tail_pad <- df[0,]
  tail_pad[nrow(tail_pad)+round(peak_window/2),] <- NA
  head_pad <- as.data.frame(head_pad)
  tail_pad <- as.data.frame(tail_pad)
  df <- rbind(head_pad,as.data.frame(df),tail_pad)
  colnames(df)[which(colnames(df) == measure1)] <- "set_measure"  
  df2 <- df %>%
    mutate(peak = rollapply(set_measure,width = peak_window,align = "center",
                            function(x) ifelse(all(is.na(x)),NA,
                                               which(x == max(x,na.rm = TRUE))==round(peak_window/2)),
                            fill = NA,partial = FALSE))
  colnames(df2)[which(colnames(df2) == "set_measure")] <- measure1 
  colnames(df2)[which(colnames(df2) == "peak")] <- paste0(measure1,"_peak_logi") 
  
  measure2 = measure[2]
  colnames(df2)[which(colnames(df2) == measure2)] <- "set_measure"  
  df3 <- df2 %>%
    mutate(peak = rollapply(set_measure,width = peak_window,align = "center",
                            function(x) ifelse(all(is.na(x)),NA,
                                               which(x == max(x,na.rm = TRUE))==round(peak_window/2)),
                            fill = NA,partial = FALSE))
  colnames(df3)[which(colnames(df3) == "set_measure")] <- measure2
  colnames(df3)[which(colnames(df3) == "peak")] <- paste0(measure2,"_peak_logi") 
  
  return(df3)
}

# df = pt_data3
# slope = c("h_slope","q_slope")
# slope_percentile = 0.95
# peak = c("h_peak_logi","q_peak_logi")
# date = "DateTime"
# test_window = 5

.pt_data4_fun <- function(df,test_window = 5,
                          slope = c("h_slope","q_slope"),
                          slope_percentile = 0.95,
                          peak = c("h_peak_logi","q_peak_logi"),
                          date = "DateTime"){
  # Set column names for analysis of first time series. 
  # Generalized to faciliate code being used for matching peaks in data other than q and h.
  colnames(df)[which(colnames(df) == slope[1])] <- "set_measure"  
  colnames(df)[which(colnames(df) == peak[1])] <- "set_peak"
  
  # Add to df a slope threshold for first series (defined by the slope_percentile argument). 
  # For those measurments marked as a peak, get slope for test_window. 
  df2 <- df %>% 
    mutate(m_thresh = quantile(abs(df$set_measure),slope_percentile,na.rm = T),
           test = ifelse(set_peak == TRUE, rollapply(set_measure,width = test_window,
                                                     function(x) ifelse(all(is.na(x)),NA,
                                                                        max(x) >= m_thresh),
                                                     fill = NA,partial = FALSE),
                         NA))
  
  # Change column names back to original name
  colnames(df2)[which(colnames(df2) == "m_thresh")] <- paste0(slope[1],slope_percentile) 
  colnames(df2)[which(colnames(df2) == "set_measure")] <- slope[1] 
  colnames(df2)[which(colnames(df2) == "set_peak")] <- peak[1]
  colnames(df2)[which(colnames(df2) == "test")] <- paste0(slope[1],"_test")
  
  # Set column names for analysis of second time series. 
  colnames(df2)[which(colnames(df2) == slope[2])] <- "set_measure"  
  colnames(df2)[which(colnames(df2) == peak[2])] <- "set_peak"  
  
  # Add to df a slope threshold for second series (defined by the slope_percentile argument). 
  # For those measurments marked as a peak, get slope for test_window. 
  df3 <- df2 %>%
    mutate(m_thresh = quantile(abs(df2$set_measure),slope_percentile,na.rm = T),
           test = ifelse(set_peak == TRUE, rollapply(set_measure,width = test_window,
                                                     function(x) ifelse(all(is.na(x)),NA,
                                                                        max(x) >= m_thresh),
                                                     fill = NA,partial = FALSE),
                         NA))
  
  # Set names back to original names.
  colnames(df3)[which(colnames(df3) == "m_thresh")] <- paste0(slope[2],slope_percentile) 
  colnames(df3)[which(colnames(df3) == "set_measure")] <- slope[2] 
  colnames(df3)[which(colnames(df3) == "set_peak")] <- peak[2]
  colnames(df3)[which(colnames(df3) == "test")] <- paste0(slope[2],"_test")
  
  return(df3)
}
# 
# df = (filter(pt_data4,date(DateTime) == "2018-03-15"))
# df = (filter(pt_data4,date(DateTime) == "2019-02-14"))
# df = filter(pt_data4,between(DateTime,ymd_hm("2019-02-14 10:00",tz=proj_tz),ymd_hm("2019-02-14 15:30",tz=proj_tz)))
# df = filter(pt_data4,between(DateTime,ymd_hm("2018-05-16 00:00",tz=proj_tz),ymd_hm("2018-05-17 00:00",tz=proj_tz)))

# df = pt_data4[1020:1070,]
# q_logi = "q_slope_test"
# h_logi = "h_slope_test"
# date = "DateTime"
# test_window = match_window
# search_direction = "right"
# x = df$q_slope_test
# df$DateTime
.match_test <- function(df,
                        test_window = 23,
                        q_logi = "q_slope_test",
                        h_logi = "h_slope_test",
                        date = "DateTime",
                        search_direction){
  colnames(df)[which(colnames(df) == q_logi)] <- "get_q"  
  colnames(df)[which(colnames(df) == h_logi)] <- "set_h"
  colnames(df)[which(colnames(df) == date)] <- "set_date"
  df2 <- df %>%
    as_tbl_time(index = set_date) %>%
    mutate(match_count = ifelse(set_h == TRUE, rollapply(get_q,align = search_direction,width = test_window,
                                                         function(x) ifelse(all(is.na(x)),NA,
                                                                            length(which(x == TRUE))),
                                                         fill = NA,partial = FALSE),
                                NA) ,
           match_row = ifelse(set_h == TRUE & match_count > 0, rollapply(get_q,align = search_direction,width = test_window,
                                                     function(x) ifelse(all(is.na(x)),NA,
                                                                        suppressWarnings(max(which(x == TRUE),na.rm = T))),
                                                     fill = NA,partial = FALSE),
                                NA))
  
  # Make match_row relative to stage peak row number, bc the rollapply reference depends on search direction:
  # for "left", match_row is relative to the stage row number.
  # for "right", match_row is relative to the match window. See below. Make both relative to stage row number.
  if(search_direction == "left"){
    df2 <- df2 %>%
      mutate(match_row = if_else(!is.na(match_row),  
                              match_row -1L,                   
                              NA_real_))                            
                                                        
    
  } else if(search_direction == "right") {
    df2 <- df2 %>%
      mutate(match_row = if_else(!is.na(match_row), 
                              (match_window-match_row),
                              NA_real_)) 
  }

  # Clean up column names
  colnames(df2)[which(colnames(df2) == "set_date")] <- date
  colnames(df2)[which(colnames(df2) == "get_q")] <- q_logi 
  colnames(df2)[which(colnames(df2) == "set_h")] <- h_logi

  
  # remove padding
  df2 <- filter(df2,!is.na(DateTime))
  
  return(df2)
}

# 
# tmp <- filter(sensor_data,between(DateTime,ymd("2017-08-05",tz=proj_tz),ymd("2017-08-06",tz=proj_tz)))
# View(select(tmp,DateTime,h_peak_logi,q_peak_logi,h_slope_test,q_slope_test,match_date23))


.override_fun <- function(dt_df,record_df) {
  for(i in 1:nrow(dt_df)){
    dt_q = dt_df[[i,1]]
    dt_h = dt_df[[i,2]]
    record_df[which(record_df$DateTime == dt_q),"q_slope_test"] <- TRUE
    record_df[which(record_df$DateTime == dt_h),"h_slope_test"] <- TRUE
  }
  return(record_df)
}

.interval_union <- function(df,int_col) {
  colnames(df)[which(colnames(df) == int_col)] <- "get_interval"  
  df2 <- df %>%
    mutate(view_interval2 = rollapply(get_interval,width = 3,align = "center",
                                      function(x){
                                        if(int_overlaps(x[2],x[1])){
                                          union(x[2],x[1])
                                        } else if(int_overlaps(x[2],x[3])){
                                          union(x[2],x[3])
                                        } else {
                                          x[2]
                                        }
                                      },
                                      fill = NA,partial = FALSE))
  return(df2)
}

###########################################################
## plotting functions for pt_data3, 4, and 5
###########################################################
.plot_data3_fun <- function(data,facet_count,window,groupscale){
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  names(cbPalette) <- as.character(expression(grey,orange,lightblue,green,yellow,darkblue,orange,purple))
  model_color <- cbPalette[c("grey","lightblue","purple","orange")]
  names(model_color) <-as.character(c("Discharge","Stage","Q Peak","H Peak"))
  
  # plot_data3 = data %>%
  #   mutate(cut_ts = cut(DateTime,facet_count),
  #          cut_ts = forcats::fct_explicit_na(cut_ts)) %>%
  #   group_by(cut_ts) %>%
  #   mutate(q_rs = scales::rescale(q),
  #          h_rs = scales::rescale(h))
  
  
  if(groupscale == TRUE){
    plot_data3 <- data %>%
      mutate(cut_ts = cut(DateTime,facet_count),
             cut_ts = forcats::fct_explicit_na(cut_ts)) %>%
      group_by(cut_ts) %>%
      mutate(q_rs = scales::rescale(q),
             h_rs = scales::rescale(h))
  } else if(groupscale == FALSE) {
    plot_data3 <- data %>%
      mutate(cut_ts = cut(DateTime,facet_count),
             cut_ts = forcats::fct_explicit_na(cut_ts)) %>%
      group_by() %>%
      mutate(q_rs = scales::rescale(q),
             h_rs = scales::rescale(h))
  }
  
  
  plot_return = ggplot(data = plot_data3) + 
    geom_line(aes(x = DateTime, y = h_rs,color = "Stage"),size = 0.5) +       
    geom_line(aes(x = DateTime, y = q_rs,color = "Discharge"),size = 0.5) +  
    geom_point(data = filter(plot_data3,h_peak_logi == TRUE),
               aes(x = DateTime, y = h_rs,color = "H Peak"),size = 1) +
    geom_point(data = filter(plot_data3,q_peak_logi == TRUE),
               aes(x = DateTime, y = q_rs,color = "Q Peak"),size = 1) +
    scale_color_manual(name = '',
                       values = model_color)+
    ylab("Rescaled h and q") +
    xlab("Date") +
    ggtitle(paste0(round(window*15/60,2)," hour peak"))+
    # slope_percentile,"% slope cutoff / ",
    # round(match_window*15/60),"hr match window")) +
    theme_minimal(base_size =12) +
    # theme_void(base_size =24) +
    theme(strip.background = element_blank(), 
          strip.text.x = element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.line = element_line(size = 0.5, linetype = "solid",
                                   colour = "black")) +
    facet_wrap(~cut_ts ,scales = "free",ncol = 1)
  return(list(plot_return,plot_data3))
}



.plot_data4_fun <- function(data,facet_count,window,groupscale=FALSE){
  plot_data4 <- filter(data,!is.na(DateTime))
  if(groupscale == TRUE){
    plot_data4 <- plot_data4 %>%
        mutate(cut_ts = cut(DateTime,facet_count),
               cut_ts = forcats::fct_explicit_na(cut_ts)) %>%
        group_by(cut_ts) %>%
        mutate(q_rs = scales::rescale(q),
               h_rs = scales::rescale(h))
  } else if(groupscale == FALSE) {
    plot_data4 <- plot_data4 %>%
      mutate(cut_ts = cut(DateTime,facet_count),
             cut_ts = forcats::fct_explicit_na(cut_ts)) %>%
      group_by() %>%
      mutate(q_rs = scales::rescale(q),
             h_rs = scales::rescale(h))
  }
  

  plot_data4 <- plot_data4 %>% 
    group_by() %>% 
    select(DateTime,cut_ts,h_rs,q_rs,h,q,h_peak_logi,q_peak_logi,h_slope_test,q_slope_test) %>% 
    mutate(q_peakmark = if_else(q_slope_test == TRUE, "q peak slope",
                                if_else(q_peak_logi == TRUE,"q peak","q")),
           h_peakmark = if_else(h_slope_test == TRUE, "h peak slope",
                                if_else(h_peak_logi == TRUE,"h peak","h")))
  
  
  hd <- plot_data4 %>% 
    select(DateTime,cut_ts,h,h_rs,h_peakmark) %>% 
    rename(meas = h,smeas = h_rs,peakmark = h_peakmark) %>% 
    mutate(meastype = "h")
  qd <- plot_data4 %>% 
    select(DateTime,cut_ts,q,q_rs,q_peakmark) %>% 
    rename(meas = q,smeas = q_rs,peakmark = q_peakmark) %>% 
    mutate(meastype = "q")
  
  plot_data4 <- bind_rows(hd,qd) %>% 
    mutate(peakmark = if_else(is.na(peakmark),meastype,peakmark))

  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#FF0000")
  names(cbPalette) <- as.character(expression(grey,orange,lightblue,green,yellow,darkblue,orange,purple,red))
  model_color <- cbPalette[c("grey","lightblue","purple","orange","orange","purple")]
  names(model_color) <-as.character(c("q","h","q peak","h peak","h peak slope","q peak slope"))
  
  model_size <- c(0,0,2,2,2,2)
  names(model_size) <-as.character(c("q","h","q peak","h peak","h peak slope","q peak slope"))
  
  model_shape <- c(NA,NA,16,16,17,17)
  names(model_shape) <-as.character(c("q","h","q peak","h peak","h peak slope","q peak slope"))
  
  plot_return <- ggplot() +
    geom_line(data = plot_data4,aes(x = DateTime, y = smeas,color = meastype),size = 0.5) +
    geom_point(data = plot_data4,aes(x = DateTime, y = smeas,
                                      color = peakmark,
                                      size = peakmark,
                                      shape = peakmark)) +
    scale_color_manual(name = '',
                       values = model_color) +
    scale_size_manual(name = '',
                      values = model_size) +
    scale_shape_manual(name = '',
                       values = model_shape) +
    
    ylab("Rescaled h and q") +
    xlab("Date") +
    ggtitle(paste0(round(window*15/60,2)," hour peak"))+
    # slope_percentile,"% slope cutoff / ",
    # round(match_window*15/60),"hr match window")) +
    theme_minimal(base_size =12) +
    # theme_void(base_size =24) +
    theme(strip.background = element_blank(), 
          strip.text.x = element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.line = element_line(size = 0.5, linetype = "solid",
                                   colour = "black")) +
    facet_wrap(~cut_ts ,scales = "free",ncol = 1)
    return(list(plot_return,plot_data4))
}

# data = pt_data5
# facet_count = 6
# window = peak_window
# groupscale = F


.plot_data5_fun <- function(data,facet_count,window,groupscale=FALSE){
  if(groupscale == TRUE){
    plot_data5 <- data %>%
      mutate(cut_ts = cut(DateTime,facet_count)) %>%
      group_by(cut_ts) %>%
      mutate(q_rs = scales::rescale(q),
             h_rs = scales::rescale(h))
  } else if(groupscale == FALSE) {
    plot_data5 <- data %>%
      mutate(cut_ts = cut(DateTime,facet_count)) %>%
      mutate(q_rs = scales::rescale(q),
             h_rs = scales::rescale(h))
  }

  plot_data5 <- plot_data5 %>%
    group_by() %>%
    select(DateTime,cut_ts,h_rs,q_rs,h,q,h_peak_logi,q_peak_logi,h_slope_test,q_slope_test,match_count,match_row,pq_row) %>% 
    mutate(h_peakmark = if_else(!is.na(match_row), "h peak match",NA_character_),
           q_peakmark = NA_character_,
           peak_n = NA_integer_)

  # num_peaks <- length(which(plot_data5$h_peakmark == "h peak match"))
  p=0
  for(i in 1:nrow(plot_data5)){
    if(!is.na(plot_data5[i,]$h_peakmark)){
      p <- p +1
      plot_data5[i,]$peak_n <- p
      plot_data5[plot_data5[i,]$pq_row,"peak_n"] <- p
      plot_data5[plot_data5[i,]$pq_row,"q_peakmark"] <- "q peak match"
    }
  }

  hd <- plot_data5 %>%
    select(DateTime,cut_ts,h,h_rs,h_peakmark,peak_n) %>%
    rename(meas = h,smeas = h_rs,peakmark = h_peakmark) %>%
    mutate(meastype = "h")
  qd <- plot_data5 %>%
    select(DateTime,cut_ts,q,q_rs,q_peakmark,peak_n) %>%
    rename(meas = q,smeas = q_rs,peakmark = q_peakmark) %>%
    mutate(meastype = "q")

  plot_data5 <- bind_rows(hd,qd) %>%
    mutate(peakmark = if_else(is.na(peakmark),meastype,peakmark))
  
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#FF0000")
  names(cbPalette) <- as.character(expression(grey,orange,lightblue,green,yellow,darkblue,orange,purple,red))
  model_color <- cbPalette[c("grey","lightblue","purple","orange","orange","purple")]
  names(model_color) <- as.character(c("q","h","h peak match","q peak match"))

  model_size <- c(0,0,2,2,2,2)
  names(model_size) <-as.character(c("q","h","h peak match","q peak match"))

  model_shape <- c(NA,NA,16,16,17,17)
  names(model_shape) <-as.character(c("q","h","h peak match","q peak match"))

  plot_return <- ggplot() +
    geom_line(data = plot_data5,aes(x = DateTime, y = smeas,color = meastype),size = 0.5) +
    geom_point(data = plot_data5,aes(x = DateTime, y = smeas,
                                     color = peakmark,
                                     size = peakmark,
                                     shape = peakmark)) +
    geom_text(data = filter(plot_data5,peakmark == "h peak match"),aes(x = DateTime,y = smeas,label = peak_n),nudge_y = 0.1) +
    scale_color_manual(name = '',
                       values = model_color) +
    scale_size_manual(name = '',
                      values = model_size) +
    scale_shape_manual(name = '',
                       values = model_shape) +

    ylab("Rescaled h and q") +
    xlab("Date") +
    ggtitle(paste0(round(window*15/60,2)," hour peak"))+
    # slope_percentile,"% slope cutoff / ",
    # round(match_window*15/60),"hr match window")) +
    theme_minimal(base_size =12) +
    # theme_void(base_size =24) +
    theme(strip.background = element_blank(),
          strip.text.x = element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.line = element_line(size = 0.5, linetype = "solid",
                                   colour = "black")) +
    facet_wrap(~cut_ts ,scales = "free",ncol = 1)
  return(list(plot_return,plot_data5))
}


.matched_peaks_fun <- function(data = pt_data6){
 
  # build dataframe of just q peaks; 
  # rename the following columns to indicate marked peak status:
  # DateTime, q, q_rs, and rnum. Note that the rnum column is changed to
  # pq_row and is the basis for the inner_join with the pt_data
  q_peaks <- data %>% 
    filter(peakmark == "q peak match") %>%
    rename(DT_pq = DateTime,pq = meas,pq_rs = smeas) %>%
    select(peak_n,DT_pq,pq,pq_rs)
  h_peaks <- data %>% 
    filter(peakmark == "h peak match") %>%
    rename(DT_ph = DateTime,ph = meas,ph_rs = smeas) %>%
    select(peak_n,DT_ph,ph,ph_rs)
  
  # Build dataframe of h peaks and q peaks, with two date columns, one for the h and one for the lag/lead q.
  matched_peaks <- bind_cols(h_peaks,q_peaks) %>% 
    # add peak number, which can be used for idetifying a particular peak
    # add time shift between the two peaks
    mutate(shift =  as.numeric(DT_ph - DT_pq,units = "hours")) %>%
    select(peak_n,DT_pq,DT_ph,shift,pq,ph,pq_rs,ph_rs)
  return(matched_peaks)
  
}

# data = pt_data6
# peak_table = matched_peaks
# view_minutes = 23

.view_port_fun <- function(data,peak_table,view_minutes){
  # Add viewing interval to matched_peaks df
  peak_table <- peak_table %>% 
    # create columns defining a period of hours before and after 
    # the marked peak to enhance QAQC visualization
    mutate(target_start = pmin(DT_ph,DT_pq), # min date between h and q date columns
           target_end = pmax(DT_ph,DT_pq), # max date between h and q date columns
           view_start = target_start - match_window*minutes(round(view_minutes/2)), # should not be more than match window
           view_end = target_end + match_window*minutes(round(view_minutes/2)),
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
  
  # Create a list of dataframes, one for each peak viewing interval. 
  # The viewing interval is defined by the view_minutes argument.
  peak_view_list <- list()
  for(i in 1:max(peak_table$peak_n)){
    peak_i <- i
    interval_i <- peak_table %>% 
      filter(peak_n == peak_i)
    interval_i <- interval_i$view_final
    view_tbl_i <- data %>% 
      filter(DateTime %within% interval_i)  
    peak_view_list[[i]] <- view_tbl_i
  }
  return(peak_view_list)
}


