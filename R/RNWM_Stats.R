RNWM_Stats<-function(dataFrame,bf_method='Eckhardt',qlimit=0,cor_method='pearson',cor_use='everything',na.rm=T){
  # This is the driver file or main Rscript
  # Make sure that the dataFrame contains model (mod) and observations (obs)

  #Assume it doesn't have missing data for now

  mod_BF<-BaseFlowSeparation(dataFrame$mod,
                                       bf_method=bf_method,
                                       k=0.93,
                                       c=quantile(streamflow,probs = 0.25, na.rm=T),
                                       filter_parameter=0.925,
                                       passes=1)
  mod_BFI<-BaseFlowIndex(mod_BF$qflow, mod_BF$bt, na.rm=T)

  mod_nev<-StormEvents(mod_BF$qft,qlimit = qlimit)

  #------------------------------------#

  obs_BF<-BaseFlowSeparation(dataFrame$obs,
                                       bf_method=bf_method,
                                       k=0.93,
                                       c=quantile(streamflow,probs = 0.25, na.rm=T),
                                       filter_parameter=0.925,
                                       passes=1)

  obs_BFI<-BaseFlowIndex(obs_BF$qflow, obs_BF$bt, na.rm=T)

  obs_nev<-StormEvents(obs_BF$qft,qlimit = qlimit)


  #-------------------------------------#
  RMSE  <-RMSE(dataFrame$mod,dataFrame$obs,na.rm=T)
  NRMSE <-PBias(dataFrame$mod,dataFrame$obs,na.rm=T)
  PBIAS <-PBias(dataFrame$mod,dataFrame$obs,na.rm=T)
  KGE   <-KGE(dataFrame$mod,dataFrame$obs, na.rm=TRUE, s.r=1, s.alpha=1, s.beta=1)
  R2    <-RSquared(dataFrame$mod,dataFrame$obs)
  NSE   <-NSE(dataFrame$mod,dataFrame$obs, nullModel=mean(dataFrame$obs, na.rm=na.rm), na.rm=TRUE)
  CORRELATION<-cor(dataFrame$mod,dataFrame$obs,method = cor_method,use=cor_use)

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
  Stat_obj$NRMSE <- NRMSE
  Stat_obj$KGE <-KGE
  Stat_obj$R2  <- R2
  Stat_obj$CORRELATION<-CORRELATION
  Stat_obj$NSE <- NSE

  return(Stat_obj)
}





