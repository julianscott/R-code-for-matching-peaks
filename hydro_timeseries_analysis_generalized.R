# hydrologic record time series analysis


# .peak_ID <- function(x){
#   return()
# }

# .sliding_mean <- function(measure,window){
#   y = tryCatch({rollapply(measure,width=window,FUN = mean,fill = NA,align = "right")},error = function(e){NA})
#   return(y)
# }

.sliding_mean <- function(measure,window){
  y = rollapply(measure,width=window,FUN = mean,fill = NA,align = "right")
  return(y)
}


# y = c(1,2,3)
# x = c(NA,NA,NA)

.rolling_lm <- rollify(.f = function(x = DateTime,y) {
  tryCatch({coef(lm(y ~ x))[2]},error = function(e){NA})
},
window = slope_window,
na_value = NA,
unlist = TRUE)

#  This is really neat way to get all of the rolling regression data 
#  But we only need the coef. 
#  See "http://lenkiefer.com/2017/10/26/predicting-recessions-with-dynamic-model-averaging/"
# tmp <- as_tbl_time(tmp,index = datetime7)
# tmp2 <- tmp %>%
#   mutate(roll_lm = rolling_lm(datetime7,q)) %>%
#   filter(!is.na(roll_lm)) %>%
#   mutate(tidied = purrr::map(roll_lm, broom::tidy)) %>%
#   unnest(tidied) %>%
#   filter(term == "x") %>%
#   select(datetime7, q, estimate) -> df4
# df4




# df <- pt_data2
# measure = c("h","q")
# peak_window = 479
# date = "DateTime"

.peak_fun <- function(df,measure = c("h","q"),peak_window = 32,date = "DateTime"){
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

.peak_test <- function(df,test_window = 5,slope = c("h_slope","q_slope"),slope_percentile = 0.95,peak = c("h_peak_logi","q_peak_logi"),date = "DateTime"){
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

# 
# q_logi = "q_slope_test"
# h_logi = "h_slope_test"
# date = "DateTime"
# test_window = match_window
# search_direction = "right"
# x = df$q_slope_test
# df$DateTime
.match_test <- function(df,test_window = 23,q_logi = "q_slope_test",h_logi = "h_slope_test",date = "DateTime",search_direction){
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
           match_row = ifelse(set_h == TRUE, rollapply(get_q,align = search_direction,width = test_window,
                                                     function(x) ifelse(all(is.na(x)),NA,
                                                                        suppressWarnings(max(which(x == TRUE),na.rm = T))),
                                                     fill = NA,partial = FALSE),
                                NA))
  
  # Make match_row relative to stage peak row number, bc the rollapply reference depends on search direction:
  # for "left", match_row is relative to the stage row number.
  # for "right", match_row is relative to the match window. See below. Make both relative to stage row number.
  if(search_dir == "left"){
    df2 <- df2 %>%
      mutate(match_row = if_else(!is.na(match_row),  
                              match_row -1L,                   
                              NA_real_))                            
                                                        
    
  } else if(search_dir == "right") {
    df2 <- df2 %>%
      mutate(match_row = if_else(!is.na(match_row), 
                              as.integer(match_window-match_row),
                              NA_integer_)) 
  }

  # Clean up column names
  colnames(df2)[which(colnames(df2) == "set_date")] <- date
  colnames(df2)[which(colnames(df2) == "get_q")] <- q_logi 
  colnames(df2)[which(colnames(df2) == "set_h")] <- h_logi
  colnames(df2)[which(colnames(df2) == "match_count")] <- paste0("match_count",test_window)
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

.plot_data3_fun <- function(data,facet_count,window){
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  names(cbPalette) <- as.character(expression(grey,orange,lightblue,green,yellow,darkblue,orange,purple))
  model_color <- cbPalette[c("grey","lightblue","purple","orange")]
  names(model_color) <-as.character(c("Discharge","Stage","Q Peak","H Peak"))
  
  plot_data3 = data %>%
    mutate(cut_ts = cut(DateTime,facet_count)) %>%
    group_by(cut_ts) %>%
    mutate(q_rs = scales::rescale(q),
           h_rs = scales::rescale(h))
  
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
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#FF0000")
  names(cbPalette) <- as.character(expression(grey,orange,lightblue,green,yellow,darkblue,orange,purple,red))
  model_color <- cbPalette[c("grey","lightblue","purple","orange","red","red")]
  names(model_color) <-as.character(c("Discharge","Stage","Q Peak","H Peak","H Peak Slope Test","Q Peak Slope Test"))
  
  if(groupscale == TRUE){
    plot_data4 = data %>%
        mutate(cut_ts = cut(DateTime,facet_count)) %>%
        group_by(cut_ts) %>%
        mutate(q_rs = scales::rescale(q),
               h_rs = scales::rescale(h))
  } else {
    plot_data4 = data %>%
      mutate(cut_ts = cut(DateTime,facet_count)) %>%
      group_by() %>%
      mutate(q_rs = scales::rescale(q),
             h_rs = scales::rescale(h))
  }
  
  plot_return = ggplot(data = plot_data4) + 
    geom_line(aes(x = DateTime, y = h_rs,color = "Stage"),size = 0.5) +       
    geom_line(aes(x = DateTime, y = q_rs,color = "Discharge"),size = 0.5) +  
    geom_point(data = filter(plot_data4,h_peak_logi == TRUE),
               aes(x = DateTime, y = h_rs,color = "H Peak"),size = 1) +
    geom_point(data = filter(plot_data4,q_peak_logi == TRUE),
               aes(x = DateTime, y = q_rs,color = "Q Peak"),size = 1) +
    geom_point(data = filter(plot_data4,h_slope_test == TRUE),
               aes(x = DateTime, y = h_rs,color = "H Peak Slope Test"),size = 3,shape = 1) +
    geom_point(data = filter(plot_data4,q_slope_test == TRUE),
               aes(x = DateTime, y = q_rs,color = "Q Peak Slope Test"),size = 3,shape = 1) +
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
  return(list(plot_return,plot_data4))
}


.plot_data5_fun <- function(data,facet_count,window,groupscale=FALSE){
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#FF0000")
  names(cbPalette) <- as.character(expression(grey,orange,lightblue,green,yellow,darkblue,orange,purple,red))
  model_color <- cbPalette[c("grey","lightblue","purple","orange","red","red")]
  names(model_color) <-as.character(c("Discharge","Stage","Q Peak","H Peak","H Peak matched","Q Peak matched"))
  
  # if(search_dir == "left"){
  #   plot_data5 <- plot_data5 %>%
  #     mutate(rnum = row_number(),
  #            sdir = search_dir,
  #            pq_row = if_else(!is.na(match_row), #    For the given row, if match_row is not NA, 
  #                             rnum + match_row -1L, # then add the row number to the number of 
  #                             NA_real_))#             rows separating the row and the row with  
  #   #             the matched peak (match_row), minus 1. Else, NA.
  #   
  # } else if(search_dir == "right") {
  #   plot_data5 <- plot_data5 %>%
  #     mutate(rnum = row_number(),
  #            sdir = search_dir,
  #            pq_row = if_else(!is.na(match_row), # Just like above, but modified for searching 'right'
  #                             rnum - as.integer(match_window-match_row),
  #                             NA_integer_)) 
  # }
  
  # rows noted as matched Q peaks
  pqrows <- plot_data5[filter(plot_data5,!is.na(pq_row))$pq_row,]
  # rows noted as matched H peaks
  phrows <- filter(plot_data5,!is.na(pq_row))
  # Bind matched peaks together, arrange by date, add column distinguising column
  mpeaks <- rbind(pqrows,phrows) %>%
    arrange(DateTime) %>% 
    mutate(PointType = if_else(is.na(match_row),"Q Peak matched","H Peak matched"))

  # Add distinguishing column to main data frame
  plot_data5$PointType <- "All Data"
  
  # combine
  tmp <- bind_rows(plot_data5,mpeaks) %>% 
    arrange(DateTime)
  
  # add plotting groups and scale
  if(groupscale == TRUE){
    tmp <- tmp %>%
      mutate(cut_ts = cut(DateTime,facet_count)) %>%
      group_by(cut_ts) %>%
      mutate(q_rs = scales::rescale(q),
             h_rs = scales::rescale(h))
  } else if(groupscale == FALSE){
    tmp <- tmp %>%
      mutate(cut_ts = cut(DateTime,facet_count)) %>%
      group_by() %>%
      mutate(q_rs = scales::rescale(q),
             h_rs = scales::rescale(h))
  }
  
  
  plot_return = ggplot() + 
    geom_line(data = filter(tmp,PointType == "All Data"),
              aes(x = DateTime, y = h_rs,color = "Stage"),size = 0.5) +       
    geom_line(data = filter(tmp,PointType == "All Data"),
              aes(x = DateTime, y = q_rs,color = "Discharge"),size = 0.5) +  
    # geom_point(data = filter(plot_data5,h_peak_logi == TRUE),
    #            aes(x = DateTime, y = h_rs,color = "H Peak"),size = 1) +
    # geom_point(data = filter(plot_data5,q_peak_logi == TRUE),
    #            aes(x = DateTime, y = q_rs,color = "Q Peak"),size = 1) +
    geom_point(data = filter(tmp,PointType == "H Peak matched"),
               aes(x = DateTime, y = h_rs,color = "H Peak matched"),size = 3,shape = 1) +
    geom_point(data = filter(tmp,PointType == "Q Peak matched"),
               aes(x = DateTime, y = q_rs,color = "Q Peak matched"),size = 3,shape = 1) +
    scale_color_manual(name = '',
                       values = model_color)+
    ylab("Rescaled h and q") +
    xlab("Date") +
    ggtitle(paste0("Matched ", round(window*15/60,2)," hour peaks"))+
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



.matched_peaks_fun <- function(df = pt_data5,
                               search_dir = search_direction,
                               match_window = match_window,
                               view_minutes = 23){
  # Depending on whether the search direction is left (sensor lags the downstream stream gage) 
  # or right (sensor leads the upstream stream gage)), this code produces three new columns:
  # 1. rnum is the row number of the datetime vector that will stay fixed to the stage record
  # 2. sdir is a confirmation indicator of the search direction that is chosen
  # 3. pq_row is the row number of a matched peak q (found within the match window, given the search direction)
  
  if(search_dir == "left"){
    df <- df %>%
      mutate(rnum = row_number(),
             sdir = search_dir,
             pq_row = if_else(!is.na(match_row), #    For the given row, if match_row is not NA, 
                              rnum + match_row -1L, # then add the row number to the number of 
                              NA_real_))#             rows separating the row and the row with  
                                        #             the matched peak (match_row), minus 1. Else, NA.
  } else {
    df <- df %>%
      mutate(rnum = row_number(),
             sdir = search_dir,
             pq_row = if_else(!is.na(match_row), # Just like above, but modified for searching 'right'
                              rnum - as.integer(match_window-match_row),
                              NA_integer_)) 
  }
  
  # scale q and h for comparing peaks and cut time series for viewing
  df <- df %>%
    mutate(q_rs = scales::rescale(q),
           h_rs = scales::rescale(h),
           cut_ts = cut(DateTime,cut_number))
  
  # Get just rows with peak discharge measurments
  pq_rows <- df[which(!is.na(df$pq_row)),]$pq_row
  
  # build dataframe of just q peaks; 
  # rename the following columns to indicate marked peak status:
  # DateTime, q, q_rs, and rnum. Note that the rnum column is changed to
  # pq_row and is the basis for the inner_join with the pt_data
  q_peaks <- select(df,rnum,DateTime,q_slope_test,q,q_rs) %>%
    filter(rnum %in% pq_rows) %>%
    rename(pq_row = rnum,DT_pq = DateTime,pq = q,pq_rs = q_rs) %>%
    select(pq_row,DT_pq,pq,pq_rs)
  
  # Build dataframe of h peaks and q peaks, with two date columns, one for the h and one for the lag/lead q.
  # this is accomplished using the inner_join with the q_peaks table based on the pq_row column.
  # Rename stage columns to indicate peak status.
  matched_peaks <- filter(df,!is.na(pq_row)) %>%
    rename(DT_ph = DateTime,ph = h,ph_rs = h_rs) %>%
    inner_join(q_peaks,by = "pq_row") %>%
    # add peak number, which can be used for idetifying a particular peak
    # add time shift between the two peaks
    mutate(peak_n = row_number(),
           cut_peak_n = cut(peak_n,2),
           shift =  as.numeric(DT_ph - DT_pq,units = "hours")) %>%
    select(peak_n,DT_pq,DT_ph,shift,pq,ph,pq_rs,ph_rs,cut_ts,cut_peak_n) %>%
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
  
  return(matched_peaks)
  
}
