BaseFlowIndex<-function(streamflow, baseflow, na.rm=T){
  BFI<-sum(baseflow,na.rm=na.rm)/sum(streamflow,na.rm=na.rm)
  return(BFI)
}
