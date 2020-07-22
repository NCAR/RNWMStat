StreamflowElasticity <- function(dfRR, WS_area){

  # dfRR: 1st column - streamflow; 2nd column - precip;
  # 3rd column: Date in POSIXct format

  tmpdf<-dfRR
  tmpdf$AreaM2<-WS_area
  tmpdf$year  <-lubridate::year(tmpdf$Date)
  tmpdf$month <-lubridate::month(tmpdf$Date)
  tmpdf$Water_Year <- NA
  tmpdf$mod_mm <-(tmpdf$mod*3600/WS_area)*1000 # to mm
  tmpdf$obs_mm <-(tmpdf$obs*3600/WS_area)*1000 # to mm
  for(i in 1:nrow(tmpdf)){
    if(tmpdf$month[i]>9){
      tmpdf$Water_Year[i]<-tmpdf$year[i]+1
    }else{
      tmpdf$Water_Year[i]<-tmpdf$year[i]
    }
  }
  gbdf <- dplyr::group_by(tmpdf,Water_Year)
  df1  <- dplyr::summarise(gbdf,sum_PCP  = sum(PCP,na.rm=T),
                           sum_qmod = sum(mod_mm,na.rm=T),
                           sum_qobs = sum(obs_mm,na.rm=T))
  df1 <- dplyr::ungroup(df1)
  dP <- diff(df1$sum_PCP)
  mP <- mean(df1$sum_PCP)
  dQm <- diff(df1$sum_qmod)
  mQm <- mean(df1$sum_qmod)
  dQo <- diff(df1$sum_qobs)
  mQo <- mean(df1$sum_qobs)

  se_mod <- quantile(dQm/dP * (mQm/mP),probs=0.5) 
  se_obs <- quantile(dQo/dP * (mQo/mP),probs=0.5) 
  l1 <- list(se_obs,se_mod)
  names(l1) <- c("se_obs","se_mod")
  return(l1)
}   

