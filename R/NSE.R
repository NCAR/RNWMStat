NSE <- function (model, observation, nullModel=mean(observation, na.rm=na.rm), na.rm=TRUE) {
  err1 <- sum((model - observation)^2, na.rm=na.rm)
  err2 <- sum((model - nullModel)^2, na.rm=na.rm)
  nse <- 1 - (err1/err2)
  return(nse)
}
