RNWM_Stats<-function(dataFrame,bf_method='Eckhardt',qlimit=0,na.rm=T){
  # This is the driver file or main Rscript
  # Make sure that the dataFrame contains model (mod) and observations (obs)

  #Assume it doesn't have missing data for now

  mod_BF<-BaseFlowSeperation(dataFrame$mod,
                                       bf_method=bf_method,
                                       k=0.93,
                                       c=quantile(streamflow,probs = 0.25, na.rm=T),
                                       filter_parameter=0.925,
                                       passes=1)
  mod_BFI<-BaseFlowIndex(mod_BF$qflow, mod_BF$bt, na.rm=T)

  mod_nev<-StormEvents(mod_BF$qft,qlimit = qlimit)

  #------------------------------------#

  obs_BF<-BaseFlowSeperation(dataFrame$obs,
                                       bf_method=bf_method,
                                       k=0.93,
                                       c=quantile(streamflow,probs = 0.25, na.rm=T),
                                       filter_parameter=0.925,
                                       passes=1)

  obs_BFI<-BaseFlowIndex(obs_BF$qflow, obs_BF$bt, na.rm=T)

  obs_nev<-StormEvents(obs_BF$qft,qlimit = qlimit)

  #-------------------------------------#
  
  dataFrame <- CalcFdc(dataFrame, colName='obs')              #Calculate flow exceedance for observed
  dataFrame <- CalcFdc(dataFrame, colName='mod')              #Calculate flow exceedance for modeled
  
  flowSplineObs <- CalcFdcSpline(dataFrame, colName = 'obs')  #Calculate the spline for observed
  flowSplineMod <- CalcFdcSpline(dataFrame, colName = 'mod')  #Calculate the spline for modeled
  
  Q10.50_Obs <- flowSplineObs(0.1)/flowSplineObs(0.5)         #Calculate the 'annual mean of the flow exceeded
  Q10.50_Mod <- flowSplineMod(0.1)/flowSplineMod(0.5)         #10% of the time' for the observed and modeled 
  
  #Calculate the slope of the center of FDC for the observed and modeled data (Boscarello paper)
  slopeFDC_Obs <- (log(flowSplineObs(0.33)) - log(flowSplineObs(0.66))) / (0.66 - 0.33)     
  slopeFDC_Mod <- (log(flowSplineMod(0.33)) - log(flowSplineMod(0.66))) / (0.66 - 0.33)
  
  #-------------------------------------#
  RMSE <-RMSE(dataFrame$mod,dataFrame$obs,na.rm=T)
  PBIAS<-PBias(dataFrame$mod,dataFrame$obs,na.rm=T)

  Stat_obj<-list()
  Stat_obj$Obs_BaseFlow<-obs_BF$bt
  Stat_obj$Obs_QuickFlow<-obs_BF$qft
  Stat_obj$Obs_BaseFlowIndex<-obs_BFI
  Stat_obj$Obs_nStormEvents<-obs_nev$nEvents

  Stat_obj$Mod_BaseFlow<-mod_BF$bt
  Stat_obj$Mod_QuickFlow<-mod_BF$qft
  Stat_obj$Mod_BaseFlowIndex<-mod_BFI
  Stat_obj$Mod_nStormEvents<-mod_nev$nEvents

  Stat_obj$RMSE  <-RMSE
  Stat_obj$PBIAS <- PBIAS
  return(Stat_obj)
}





