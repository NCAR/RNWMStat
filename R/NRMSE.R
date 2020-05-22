NRMSE <- function (model, observation, na.rm=TRUE) {
  err <- sum((model - observation)^2, na.rm=na.rm)/(min(sum(!is.na(model)),sum(!is.na(observation))))
  rmserr <- sqrt(err) / ( max(observation, na.rm=na.rm) - min(observation, na.rm=na.rm) ) * 100
  return(rmserr)
}
