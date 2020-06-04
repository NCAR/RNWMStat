RNWM_Stats<-function(dataFrame,WS_area,bf_method='Lyne-Nathan',qlimit=0,cor_use='everything',na.rm=T){

  # This is the driver file or main Rscript
  # Make sure that the dataFrame contains date, model (mod) and observations (obs)

  #Assume it doesn't have missing data for now

  dataFrame<-na.omit(dataFrame)

  #------------------------------------#
  # Model calculations
  #------------------------------------#

    #------------------------------------#
    # Calculate baseflow (mod_BF)
    #------------------------------------#

    mod_BF<-BaseFlowSeparation(dataFrame$mod,
                                       bf_method=bf_method,
                                       k=0.93,
                                       c=quantile(streamflow,probs = 0.25, na.rm=T),
                                       filter_parameter=0.925,
                                       passes=1)

    #------------------------------------#
    # Calculate baseflow index (mod_BFI)
    #------------------------------------#

    mod_BFI<-BaseFlowIndex(mod_BF$qflow, mod_BF$bt, na.rm=T)

    #------------------------------------#
    # Calculate events (mod_nev)
    #------------------------------------#

    mod_nev<-StormEvents(mod_BF$qft,qlimit = qlimit)

    #------------------------------------#
    # Flow duration curve (FDC)
    #------------------------------------#

    mod_Fdc        <- CalcFdc(dataFrame, colName='mod')              #Calculate flow exceedance for modeled
    mod_flowSpline <- CalcFdcSpline(mod_Fdc, colName = 'mod')      #Calculate the spline for modeled
    mod_Q10.50     <- Q10_50(mod_flowSpline)                        #10% of the time' for the observed and modeled

    #Calculate the slope of the center of FDC for the observed and modeled data (Boscarello paper)
    mod_slopeFDC <- Calc_FDC_Slope(mod_flowSpline)

    #------------------------------------#
    # Event Identification
    #------------------------------------#

    Date       <-as.POSIXct(dataFrame$Date)
    mod        <-dataFrame$mod
    mod_dd     <-as.data.frame(Date)
    mod_dd$mod <-mod

    mod_dd<-na.omit(mod_dd)
    mod_evdf<-eventIdentification(mod_dd,tag='aaa',nwinSpan=36,plot=FALSE,
                            probPeak=0.9,probLowflow=0.5,probRange=0.3,
                            minEventDuration=6,maxRiseDuration=120,maxRecessionDuration=240,
                            nwinShift=12,minInterval=6,minLengthData=6)

    mod_total_events     <- nrow(mod_evdf)

    #------------------------------------#
    # Lag Time
    #------------------------------------#
    mod_LagTime<-LagTime(dataFrame$PCP,dataFrame$mod,lag.max=7)

  #------------------------------------#
  # Obs calculations
  #------------------------------------#

    #------------------------------------#
    # Calculate baseflow (obs_BF)
    #------------------------------------#

    obs_BF<-BaseFlowSeparation(dataFrame$obs,
                                       bf_method=bf_method,
                                       k=0.93,
                                       c=quantile(streamflow,probs = 0.25, na.rm=T),
                                       filter_parameter=0.925,
                                       passes=1)

    #------------------------------------#
    # Calculate baseflow index (obs_BFI)
    #------------------------------------#

    obs_BFI<-BaseFlowIndex(obs_BF$qflow, obs_BF$bt, na.rm=T)

    #------------------------------------#
    # Calculate events (obs_nev)
    #------------------------------------#

    obs_nev<-StormEvents(obs_BF$qft,qlimit = qlimit)

    #------------------------------------#
    # Flow duration curve (FDC)
    #------------------------------------#

    obs_Fdc        <- CalcFdc(dataFrame, colName='obs')              #Calculate flow exceedance for observed
    obs_flowSpline <- CalcFdcSpline(obs_Fdc, colName = 'obs')      #Calculate the spline for observed
    obs_Q10.50     <- Q10_50(obs_flowSpline)                        #Calculate the 'annual mean of the flow exceeded

    #Calculate the slope of the center of FDC for the observed and modeled data (Boscarello paper)
    obs_slopeFDC <- Calc_FDC_Slope(obs_flowSpline)

    #------------------------------------#
    # Event Identification
    #------------------------------------#

    Date       <-as.POSIXct(dataFrame$Date)
    obs        <-dataFrame$obs
    obs_dd     <-as.data.frame(Date)
    obs_dd$Obs <-obs

    obs_dd<-na.omit(obs_dd)
    obs_evdf<-eventIdentification(obs_dd,tag='aaa',nwinSpan=36,plot=FALSE,
                                  probPeak=0.9,probLowflow=0.5,probRange=0.3,
                                  minEventDuration=6,maxRiseDuration=120,maxRecessionDuration=240,
                                  nwinShift=12,minInterval=6,minLengthData=6)

    obs_total_events     <- nrow(obs_evdf)

    #------------------------------------#
    # Event Matching
    #------------------------------------#

    event_stats <- match_events_obs_mod(obs_evdf,mod_evdf,dataFrame)

    #------------------------------------#
    # Lag Time
    #------------------------------------#
    obs_LagTime<-LagTime(dataFrame$PCP,dataFrame$obs,lag.max=7)

    #------------------------------------#
    # Lag Time
    #------------------------------------#

    ROffRatio <- RunOffRatio(dataFrame,WS_area)

  #------------------------------------#
  # Metrics requiring model and obs
  #------------------------------------#

    #------------------------------------#
    # Calculate root mean square error (RMSE)
    #------------------------------------#

    RMSE <-RMSE(dataFrame$mod,dataFrame$obs,na.rm=T)

    #------------------------------------#
    # Calculate normalized root mean square error (NRMSE)
    #------------------------------------#

    NRMSE <-NRMSE(dataFrame$mod,dataFrame$obs,na.rm=T)

    #------------------------------------#
    # Calculate percent bias (PBIAS)
    #------------------------------------#

    PBIAS<-PBias(dataFrame$mod,dataFrame$obs,na.rm=T)

    #------------------------------------#
    # Calculate NSE (NSE)
    #------------------------------------#
    NSE   <-NSE(dataFrame$mod,dataFrame$obs, nullModel=mean(dataFrame$obs, na.rm=na.rm), na.rm=TRUE)


    #------------------------------------#
    # Calculate KGE (KGE)
    #------------------------------------#
    KGE   <-KGE(dataFrame$mod,dataFrame$obs, na.rm=TRUE, s.r=1, s.alpha=1, s.beta=1)


    #------------------------------------#
    # Calculate correlation (R2)
    #------------------------------------#
    R2    <-RSquared(dataFrame$mod,dataFrame$obs)


    #------------------------------------#
    # Calculate Pearson R (Pearson_R)
    #------------------------------------#
    COR_Pearson<-cor(dataFrame$mod,dataFrame$obs,method ='pearson',use=cor_use)

    #------------------------------------#
    # Calculate Spearman R (Spearman_R)
    #------------------------------------#
    COR_Spearman<-cor(dataFrame$mod,dataFrame$obs,method ='spearman',use=cor_use)

    #------------------------------------#
    # Calculate Spearman R (Spearman_R)
    #------------------------------------#

  #------------------------------------#
  # Populate statistics object
  #------------------------------------#

  Stat_obj<-list()
  Stat_obj$Obs_BaseFlow     <-obs_BF$bt
  Stat_obj$Obs_QuickFlow    <-obs_BF$qft
  Stat_obj$Obs_BaseFlowIndex<-obs_BFI
  Stat_obj$Obs_nStormEvents <-obs_nev$nEvents
  Stat_obj$Obs_Fdc          <-obs_Fdc
  Stat_obj$Obs_flowSpline   <-obs_flowSpline
  Stat_obj$Obs_Q10.50       <-obs_Q10.50
  Stat_obj$Obs_slopeFDC     <-obs_slopeFDC

  Stat_obj$Mod_BaseFlow     <-mod_BF$bt
  Stat_obj$Mod_QuickFlow    <-mod_BF$qft
  Stat_obj$Mod_BaseFlowIndex<-mod_BFI
  Stat_obj$Mod_nStormEvents <-mod_nev$nEvents
  Stat_obj$Mod_Fdc          <-obs_Fdc
  Stat_obj$Mod_flowSpline   <-obs_flowSpline
  Stat_obj$Mod_Q10.50       <-obs_Q10.50
  Stat_obj$Mod_slopeFDC     <-obs_slopeFDC

  Stat_obj$RMSE  <-RMSE
  Stat_obj$PBIAS <-PBIAS
  Stat_obj$NRMSE <-NRMSE
  Stat_obj$KGE   <-KGE
  Stat_obj$R2    <-R2
  Stat_obj$COR_Pearson <-COR_Pearson
  Stat_obj$COR_Spearman<-COR_Spearman
  Stat_obj$NSE   <-NSE

  Stat_obj$event_stats      <-event_stats
  Stat_obj$obs_total_events <-obs_total_events
  Stat_obj$mod_total_events <-mod_total_events
  Stat_obj$obs_LagTime      <-obs_LagTime
  Stat_obj$mod_LagTime      <-mod_LagTime
  Stat_obj$ROffRatio        <-ROffRatio

  return(Stat_obj)
}





