match_events_obs_mod <- function(obs_evdf,mod_evdf,dataFrame){

  mod_event_start_hour <- as.numeric(mod_evdf$start)/(60*60)
  mod_event_peak_hour  <- as.numeric(mod_evdf$peak)/(60*60)
  mod_event_end_hour   <- as.numeric(mod_evdf$end)/(60*60)
  mod_event_length     <- mod_event_end_hour - mod_event_start_hour
  mod_total_events     <- nrow(mod_evdf)

  obs_event_start_hour <- as.numeric(obs_evdf$start)/(60*60)
  obs_event_peak_hour  <- as.numeric(obs_evdf$peak)/(60*60)
  obs_event_end_hour   <- as.numeric(obs_evdf$end)/(60*60)
  obs_event_length     <- obs_event_end_hour - obs_event_start_hour
  obs_total_events     <- nrow(obs_evdf)


  obs_evdf$start<-lubridate::force_tz(as.POSIXct(obs_evdf$start,format="%Y-%m-%d %H:%M:%OS"),tz="MST")
  obs_evdf$peak <-lubridate::force_tz(as.POSIXct(obs_evdf$peak,format="%Y-%m-%d %H:%M:%OS"),tz="MST")
  obs_evdf$end  <-lubridate::force_tz(as.POSIXct(obs_evdf$end,format="%Y-%m-%d %H:%M:%OS"),tz="MST")

  mod_evdf$start<-lubridate::force_tz(as.POSIXct(mod_evdf$start,format="%Y-%m-%d %H:%M:%OS"),tz="MST")
  mod_evdf$peak <-lubridate::force_tz(as.POSIXct(mod_evdf$peak,format="%Y-%m-%d %H:%M:%OS"),tz="MST")
  mod_evdf$end  <-lubridate::force_tz(as.POSIXct(mod_evdf$end,format="%Y-%m-%d %H:%M:%OS"),tz="MST")

  dataFrame$Date<-lubridate::force_tz(as.POSIXct(dataFrame$Date),tz="MST")

  obs_qpeak    <-vector(length = nrow(obs_evdf))
  obs_peak_time<-vector(length = nrow(obs_evdf))
  obs_ev_dur   <-vector(length = nrow(obs_evdf))

  mod_qpeak    <-vector(length = nrow(obs_evdf))
  mod_peak_time<-vector(length = nrow(obs_evdf))
  mod_ev_dur   <-vector(length = nrow(obs_evdf))

  PeakTimeDiff <-vector(length = nrow(obs_evdf))

  i<-10
  j<-10

  for (i in 1:nrow(obs_evdf)){

    for (j in 1:nrow(mod_evdf)){

      if( (mod_event_peak_hour[j]>obs_event_peak_hour[i]-6) && (mod_event_peak_hour[j]<obs_event_peak_hour[i]+6)  ){

        # print(paste(i,j,obs_event_peak_hour[i]-6,
        #             mod_event_peak_hour[j],
        #             obs_event_peak_hour[i]+6,
        #             dataFrame$obs[which(dataFrame$Date==obs_evdf$peak[i])]))

        obs_qpeak[i]     <-dataFrame$obs[which(dataFrame$Date==obs_evdf$peak[i])]
        obs_peak_time[i] <-as.numeric(obs_evdf$peak[i])
        obs_ev_dur[i]    <-(as.numeric(obs_evdf$end[i])-as.numeric(obs_evdf$start[i]))/(3600)


        mod_qpeak[i]     <-dataFrame$mod[which(dataFrame$Date==mod_evdf$peak[j])]
        mod_peak_time[i] <-as.numeric(mod_evdf$peak[j])
        mod_ev_dur[i]    <-(as.numeric(mod_evdf$end[j])-as.numeric(mod_evdf$start[j]))/(3600)


        PeakTimeDiff[i] <-  (mod_peak_time[i]-obs_peak_time[i])/3600

      }


    }

  }

  Events_match<-data.frame(obs_qpeak,mod_qpeak,
                           mod_peak_time,obs_peak_time,
                           mod_ev_dur,obs_ev_dur,
                           PeakTimeDiff)


  Events_match$obs_qpeak[Events_match$obs_qpeak==0.]<-NA

  Events_match<-na.omit(Events_match)

  return(Events_match)
}
