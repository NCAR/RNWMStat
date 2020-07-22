RNWM_Stats_nonEvent<-function(dataFrame,WS_area,bf_method='Lyne-Nathan',qlimit=0,cor_use='everything',na.rm=T){

  # Script that computes all metrics that are not event-based and collects all stats
  # into a list.  
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
                                       c=quantile(dataFrame$mod,probs = 0.25, na.rm=T),
                                       filter_parameter=0.925,
                                       passes=1)

    #------------------------------------#
    # Calculate baseflow index (mod_BFI)
    #------------------------------------#

    mod_BFI<-BaseFlowIndex(mod_BF$qflow, mod_BF$bt, na.rm=T)


    #------------------------------------#
    # Flow duration curve (FDC)
    #------------------------------------#

    mod_Fdc        <- CalcFdc(dataFrame, colName='mod')              #Calculate flow exceedance for modeled
    mod_flowSpline <- CalcFdcSpline(mod_Fdc, colName = 'mod')      #Calculate the spline for modeled
    mod_Q10.50     <- Q10_50(mod_flowSpline)                        #10% of the time' for the observed and modeled

    #Calculate the slope of the center of FDC for the observed and modeled data (Boscarello paper)
    mod_slopeFDC <- Calc_FDC_Slope(mod_flowSpline)


    #------------------------------------#
    # Lag Time
    #------------------------------------#
    mod_LagTime<-calcLagTime(dataFrame$PCP,dataFrame$mod)

    #------------------------------------#
    # (mean) runoff coefficient 
    #------------------------------------#
    mod_RunoffCoeff <- RunoffCoefficient(dataFrame[,c("mod","PCP")],WS_area)

    #------------------------------------#
    # runoff ratio by water year for model and obs 
    #------------------------------------#
    tmp <- RunOffRatio(dataFrame,WS_area)
    mod_RunoffRatioWY <- tmp$mod_runoff_ratio
    obs_RunoffRatioWY <- tmp$obs_runoff_ratio

    #------------------------------------#
    # streamflow elasticity for model and obs
    #------------------------------------#
    tmp <- StreamflowElasticity(dataFrame,WS_area)
    mod_flowElast <- tmp$se_mod
    obs_flowElast <- tmp$se_obs

  #------------------------------------#
  # Obs calculations
  #------------------------------------#

    #------------------------------------#
    # Calculate baseflow (obs_BF)
    #------------------------------------#

    obs_BF<-BaseFlowSeparation(dataFrame$obs,
                                       bf_method=bf_method,
                                       k=0.93,
                                       c=quantile(dataFrame$obs,probs = 0.25, na.rm=T),
                                       filter_parameter=0.925,
                                       passes=1)

    #------------------------------------#
    # Calculate baseflow index (obs_BFI)
    #------------------------------------#

    obs_BFI<-BaseFlowIndex(obs_BF$qflow, obs_BF$bt, na.rm=T)

    #------------------------------------#
    # Flow duration curve (FDC)
    #------------------------------------#

    obs_Fdc        <- CalcFdc(dataFrame, colName='obs')              #Calculate flow exceedance for observed
    obs_flowSpline <- CalcFdcSpline(obs_Fdc, colName = 'obs')      #Calculate the spline for observed
    obs_Q10.50     <- Q10_50(obs_flowSpline)                        #Calculate the 'annual mean of the flow exceeded

    #Calculate the slope of the center of FDC for the observed and modeled data (Boscarello paper)
    obs_slopeFDC <- Calc_FDC_Slope(obs_flowSpline)

    #------------------------------------#
    # Lag Time
    #------------------------------------#
    obs_LagTime<-calcLagTime(dataFrame$PCP,dataFrame$obs)

    #------------------------------------#
    # (mean) runoff coefficient 
    #------------------------------------#
    obs_RunoffCoeff <- RunoffCoefficient(dataFrame[,c("obs","PCP")],WS_area)

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
    # Calculate NSE of log-tansformed flow (NseLog)
    #------------------------------------#
    NseLog <- NseLog(dataFrame$mod,dataFrame$obs)

    #------------------------------------#
    # Calculate weighted NSE: w1*NSE + (1-w1)*NseLog 
    #------------------------------------#
    NseWt <- NseWt(dataFrame$mod,dataFrame$obs)

    #------------------------------------#
    # Calculate Multi-scale objective function (Msof)
    #------------------------------------#
    Msof <- Msof(dataFrame$mod,dataFrame$obs)

    #------------------------------------#
    # Calculate hyper resolution objective function (hyperResMultiObj)
    #------------------------------------#
    hyperResMultiObj <- hyperResMultiObj(dataFrame$mod,dataFrame$obs)


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
    # Calculate r1 (lognormal esitmator of correlation based on Stedinger's (1981)
    #------------------------------------#
    df1 <- noZeroFunction(dataFrame$mod,dataFrame$obs)
    mod1 <- df1[[1]]; obs1 <- df1[[2]]
    COR_r1 <- r1(mod1, obs1)

    #------------------------------------#
    # Calculate LBEm functions for NSE and KGE
    #------------------------------------#
    LBEms <- LBEms_function(mod1, obs1,period=lubridate::month(dataFrame$Date),
        monsample=720)
    LBEm_NSE <- LBEms[1]
    LBEm_KGE <- LBEms[2]

  #------------------------------------#
  # Populate statistics object
  #------------------------------------#

  Stat_obj<-list()
  Stat_obj$Obs_BaseFlow     <-obs_BF$bt
  Stat_obj$Obs_QuickFlow    <-obs_BF$qft
  Stat_obj$Obs_BaseFlowIndex<-obs_BFI
  Stat_obj$Obs_Fdc          <-obs_Fdc
  Stat_obj$Obs_flowSpline   <-obs_flowSpline
  Stat_obj$Obs_Q10.50       <-obs_Q10.50
  Stat_obj$Obs_slopeFDC     <-obs_slopeFDC
  Stat_obj$Obs_RunoffCoeff  <-obs_RunoffCoeff 
  Stat_obj$Obs_RunoffRatioWY <-obs_RunoffRatioWY 
  Stat_obj$Obs_LagTime      <-obs_LagTime
  Stat_obj$Obs_flowElast    <-obs_flowElast

  Stat_obj$Mod_BaseFlow     <-mod_BF$bt
  Stat_obj$Mod_QuickFlow    <-mod_BF$qft
  Stat_obj$Mod_BaseFlowIndex<-mod_BFI
  Stat_obj$Mod_Fdc          <-mod_Fdc
  Stat_obj$Mod_flowSpline   <-mod_flowSpline
  Stat_obj$Mod_Q10.50       <-mod_Q10.50
  Stat_obj$Mod_slopeFDC     <-mod_slopeFDC
  Stat_obj$Mod_RunoffCoeff  <-mod_RunoffCoeff 
  Stat_obj$Mod_RunoffRatioWY <-mod_RunoffRatioWY 
  Stat_obj$Mod_LagTime      <-mod_LagTime
  Stat_obj$Mod_flowElast    <-mod_flowElast

  Stat_obj$RMSE  <-RMSE
  Stat_obj$PBIAS <-PBIAS
  Stat_obj$NRMSE <-NRMSE
  Stat_obj$KGE   <-KGE
  Stat_obj$R2    <-R2
  Stat_obj$COR_Pearson <-COR_Pearson
  Stat_obj$COR_Spearman<-COR_Spearman
  Stat_obj$COR_r1 <- COR_r1
  Stat_obj$NSE   <-NSE
  Stat_obj$NseLog <- NseLog
  Stat_obj$NseWt <- NseWt
  Stat_obj$LBEm_NSE <- LBEm_NSE
  Stat_obj$LBEm_KGE <- LBEm_KGE
  Stat_obj$Msof <- Msof
  Stat_obj$hyperResMultiObj <- hyperResMultiObj
  

  return(Stat_obj)
}





