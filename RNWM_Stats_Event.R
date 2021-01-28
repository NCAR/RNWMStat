RNWM_Stats_Event<-function(dataFrame,WS_area,snow=FALSE, slow=FALSE, threshold,lagTime=24, na.rm=T) {

  # Computes all event-based metrics and collects stats into a list.
  # Make sure that the dataFrame contains date (Date, POSIXct format), model (mod),
  #  observations (obs), and precipitation amount (PCP)

  #Assume it doesn't have missing data for now
  #dataFrame<-na.omit(dataFrame)

  # identify events for observed streamflow
  listObs <- eventIdentification(dataFrame[,c("Date","obs")],snow,slow,threshold)
  data_obs <- listObs[["dataAll"]]
  events_obs0 <- listObs[["eventsAll"]]
      
  # identify events for model streamflow
  listMod <- eventIdentification(dataFrame[,c("Date","mod")],snow,slow,threshold)
  data_mod <- listMod[["dataAll"]]
  events_mod0 <- listMod[["eventsAll"]]  

  # match observed events with model events
  matchResults <- matchEvents(data_obs,data_mod,events_obs0,events_mod0,threshold,snow)
  events_obs <- matchResults[["events_obs_match"]]
  events_mod <- matchResults[["events_mod_match"]]
  ne1 <- nrow(events_obs)

  ##############  event metrics ################

  Stat_obj <- list()
  Stat_obj$matchResults <- matchResults
  if (ne1>0) {

  peak_obs <- data_obs$value[match(events_obs$peak,data_obs$time)]
  peak_mod <- data_mod$value[match(events_mod$peak,data_mod$time)]

  Stat_obj$nHit <- sum(peak_mod>=threshold)
  Stat_obj$nMiss <- sum(peak_mod<threshold)
  Stat_obj$hitRate <- Stat_obj$nHit/ne1*100
  
  # false alarms
  ix1 <- data_mod$value[match(events_mod0$peak, data_mod$time)]>=threshold
  events_mod0 <- events_mod0[ix1,]
  dates1 <- NULL
  for (i1 in 1:ne1) dates1 <- c(dates1,seq(events_mod$start[i1],events_mod$end[i1],by="hour"))
  ix1 <- which(! events_mod0$peak %in% dates1)
  Stat_obj$nFalseAlarm <- length(ix1)
  Stat_obj$falseAlarmRate <- Stat_obj$nFalseAlarm/ne1*100

  Stat_obj$timingError <- as.integer(difftime(events_mod$peak,events_obs$peak,units="hour"))
  Stat_obj$pBiasPeak <- (peak_mod-peak_obs)/peak_obs*100

  volume_mod <- sapply(1:ne1, function(x) sum(subset(data_mod,time %in% seq(events_mod$start[x],events_mod$end[x],by="hour"))$value))
  volume_obs <- sapply(1:ne1, function(x) sum(subset(data_obs,time %in% seq(events_obs$start[x],events_obs$end[x],by="hour"))$value))
  Stat_obj$pBiasVolume <- (volume_mod-volume_obs)/volume_obs*100

  # total (event+antecedent) precipitation for computing runoff ratio
  pcp_mod <- sapply(1:ne1,function(x) sum(subset(dataFrame,Date %in% seq(events_mod$start[x]-lagTime*3600,events_mod$end[x],by="hour"))$PCP,na.rm=T))
  pcp_obs <- sapply(1:ne1,function(x) sum(subset(dataFrame,Date %in% seq(events_obs$start[x]-lagTime*3600,events_obs$end[x],by="hour"))$PCP,na.rm=T))
  Stat_obj$Mod_eventRR <- volume_mod*3600/WS_area*1000/pcp_mod
  Stat_obj$Obs_eventRR <- volume_obs*3600/WS_area*1000/pcp_obs


  Stat_obj$meanTimeErr <- mean(abs(Stat_obj$timingError))
  Stat_obj$meanPeakBias <- mean(abs(Stat_obj$pBiasPeak))
  Stat_obj$meanVolumeBias <- mean(abs(Stat_obj$pBiasVolume))
  Stat_obj$Mod_meanEventRR <- mean(abs(Stat_obj$Mod_eventRR))
  Stat_obj$Obs_meanEventRR <- mean(abs(Stat_obj$Obs_eventRR))
  }

  return(Stat_obj)
}

