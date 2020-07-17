RNWM_Stats_Event<-function(dataFrame,WS_area,snow=FALSE, groundwater=FALSE, slow=FALSE, lagTime=24, na.rm=T) {

  # Script that computes all metrics that are not event-based and collects all stats
  # into a list.
  # Make sure that the dataFrame contains date (Date, POSIXct format), model (mod),
  #  observations (obs), and precipitation amount (PCP)

  #Assume it doesn't have missing data for now
  dataFrame<-na.omit(dataFrame)

  # set up parameters  
  nwinSpan1 <- 36
  minEventDuration1 <- 12
  maxRiseDuration1 <- 2*24
  maxRecessionDuration1 <- 5*24
  maxPeakDistCompound1 <- 5*24
  maxDist1 <- 2*24

  # adjust for seasonal snow
  if (snow) {
    nwinSpan1 <- 240
    minEventDuration1 <- 36
    maxRiseDuration1 <- 30*24
    maxRecessionDuration1 <- 60*24
    maxPeakDistCompound1 <- 15*24
    maxDist1 <- 7*24
  }

  # adjust for groundwater
  if (groundwater) {
    nwinSpan1 <- 72
    maxRiseDuration1 <- 30*24
    maxRecessionDuration1 <- 30*24
  }

  # adjust for slow rising and recession limbs 
  if (slow) {
    maxRiseDuration1 <- 5*24
    maxRecessionDuration1 <- 8*24
  }

  # identify events for observed streamflow
  listObs <- eventIdentification(dataFrame[,c("Date","obs")],
      nwinSpan=nwinSpan1, 
      minEventDuration=minEventDuration1, 
      maxRiseDuration=maxRiseDuration1, 
      maxRecessionDuration=maxRecessionDuration1,
      maxPeakDistCompound=maxPeakDistCompound1)

  # identify events for model streamflow
  listMod <- eventIdentification(dataFrame[,c("Date","mod")],
      nwinSpan=nwinSpan1,
      minEventDuration=minEventDuration1,
      maxRiseDuration=maxRiseDuration1,
      maxRecessionDuration=maxRecessionDuration1,
      maxPeakDistCompound=maxPeakDistCompound1)

  # match observed events with model events
  eventsMatched <- matchEvents(dataFrame[,c("Date","mod")],
     listMod[[1]], listMod[[2]], listObs[[2]], maxDist=maxDist1)

  ##############  event metrics ################

  #data from eventIdentification with short gaps filled
  obsDf <- listObs[[3]] 
  modDf <- listMod[[3]]

  Stat_obj <- list()
  Stat_obj$Obs_nEvents <- nrow(listObs[[2]])
  Stat_obj$Mod_nEvents <- nrow(listMod[[2]])
  Stat_obj$Obs_peak <- obsDf$value[match(eventsMatched$peak_obs, obsDf$time)]
  Stat_obj$Mod_peak <- modDf$value[match(eventsMatched$peak_mod, modDf$time)]

  if (nrow(eventsMatched)==0) {
    
  Stat_obj$nHit <- 0
  Stat_obj$nMiss <- Stat_obj$Obs_nEvents
  Stat_obj$nFalseAlarm <- Stat_obj$Mod_nEvents

  } else {
  
  Stat_obj$nHit <- sum(eventsMatched$match %in% c(1,2))
  Stat_obj$nMiss <- sum(eventsMatched$match==3)
  Stat_obj$nFalseAlarm <- Stat_obj$Mod_nEvents - sum(eventsMatched$match==1)
  Stat_obj$timingError <- as.integer(difftime(eventsMatched$peak_mod,eventsMatched$peak_obs,units="hour"))
  Stat_obj$pBiasPeak <- (Stat_obj$Obs_peak - Stat_obj$Mod_peak)/Stat_obj$Obs_peak*100

  dfMetric <- plyr::ldply(1:nrow(eventsMatched), function(x) {
    k1 <- match(eventsMatched$start_obs[x], obsDf$time)
    k2 <- match(eventsMatched$end_obs[x], obsDf$time)
    obs1 <- obsDf$value[k1:k2]
    k1 <- match(eventsMatched$start_mod[x], modDf$time)
    k2 <- match(eventsMatched$end_mod[x], modDf$time)
    mod1 <- modDf$value[k1:k2]

    volume_bias <- (sum(mod1)-sum(obs1))/sum(obs1)*100 

    # total (event+antecedent) precipitation for computing runoff ratio
    df0 <- subset(dataFrame,Date %in% seq(eventsMatched$start_obs[x]-lagTime*3600,eventsMatched$end_obs[x],by="hour"))
    pcp0 <- sum(df0$PCP,na.rm=T)
    
    # event runoff ratio
    mod_eventRR <- sum(mod1*3600/WS_area*1000)/pcp0
    obs_eventRR <- sum(obs1*3600/WS_area*1000)/pcp0

    c(volume_bias,mod_eventRR,obs_eventRR)})

  Stat_obj$pBiasVolume <- dfMetric[[1]]
  Stat_obj$Mod_eventRR <- dfMetric[[2]]
  Stat_obj$Obs_eventRR <- dfMetric[[3]]

  Stat_obj$meanTimeErr <- mean(Stat_obj$timingError)
  Stat_obj$meanPeakBias <- mean(Stat_obj$pBiasPeak)
  Stat_obj$meanVolumeBias <- mean(Stat_obj$pBiasVolume)
  Stat_obj$Mod_meanEventRR <- mean(Stat_obj$Mod_eventRR)
  Stat_obj$Obs_meanEventRR <- mean(Stat_obj$Obs_eventRR)
  }

  return(Stat_obj)
}

