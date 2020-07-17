RunoffCoefficient<-function(dfRR,WS_area){
  # dfRR: 1st column - streamflow; 2nd column - precip
  names(dfRR) <- c("Q","PCP")
  dfRR$Q <-(dfRR$Q*3600/WS_area)*1000 # to mm
  rc <-  sum(dfRR$Q)/sum(dfRR$PCP)
  rc
}

