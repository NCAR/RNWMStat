RunOffRatio<-function(dfRR,WS_area){
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

  sb   <- dplyr::summarise(gbdf,sum_PCP  = sum(PCP,na.rm=T),
                           sum_qmod = sum(mod_mm,na.rm=T),
                           sum_qobs = sum(obs_mm,na.rm=T),
                           obs_runoff_ratio = sum(obs_mm,na.rm=T)/sum(PCP,na.rm=T),
                           mod_runoff_ratio = sum(mod_mm,na.rm=T)/sum(PCP,na.rm=T))
  SFdb<-dplyr::ungroup(sb)
  return(SFdb)
}

