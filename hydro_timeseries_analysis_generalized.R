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
  colnames(df)[which(colnames(df) == slope[1])] <- "set_measure"  
  colnames(df)[which(colnames(df) == peak[1])] <- "set_peak"  
  df2 <- df %>%
    mutate(m_thresh = quantile(abs(df$set_measure),slope_percentile,na.rm = T),
           test = ifelse(set_peak == TRUE, rollapply(set_measure,width = test_window,
                                                     function(x) ifelse(all(is.na(x)),NA,
                                                                        max(x) >= m_thresh),
                                                     fill = NA,partial = FALSE),
                         NA))
  colnames(df2)[which(colnames(df2) == "m_thresh")] <- paste0(slope[1],slope_percentile) 
  colnames(df2)[which(colnames(df2) == "set_measure")] <- slope[1] 
  colnames(df2)[which(colnames(df2) == "set_peak")] <- peak[1]
  colnames(df2)[which(colnames(df2) == "test")] <- paste0(slope[1],"_test")
  
  colnames(df2)[which(colnames(df2) == slope[2])] <- "set_measure"  
  colnames(df2)[which(colnames(df2) == peak[2])] <- "set_peak"  
  df3 <- df2 %>%
    mutate(m_thresh = quantile(abs(df2$set_measure),slope_percentile,na.rm = T),
           test = ifelse(set_peak == TRUE, rollapply(set_measure,width = test_window,
                                                     function(x) ifelse(all(is.na(x)),NA,
                                                                        max(x) >= m_thresh),
                                                     fill = NA,partial = FALSE),
                         NA))
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
  colnames(df2)[which(colnames(df2) == "set_date")] <- date
  colnames(df2)[which(colnames(df2) == "get_q")] <- q_logi 
  colnames(df2)[which(colnames(df2) == "set_h")] <- h_logi
  colnames(df2)[which(colnames(df2) == "match_count")] <- paste0("match_count",test_window)
  colnames(df2)[which(colnames(df2) == "match_row")] <- paste0("match_row",test_window)
  # tmp <- filter(df2,between(DateTime,ymd("2017-08-05",tz=proj_tz),ymd("2017-08-07",tz=proj_tz)))
  # View(select(tmp,DateTime,h_peak_logi,q_peak_logi,h_slope_test,q_slope_test,match_date23))
  
  # df3 <- select(df2,-match_test)
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

  


