# Percent Bias
PBias<- function (model, observation, na.rm=TRUE) {
  pbs<-sum(model - observation, na.rm=na.rm)/sum(observation, na.rm=na.rm) * 100
  return(pbs)
}
