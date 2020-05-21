RMSE <- function (model, observation, na.rm=TRUE) {
  RMSError <- sqrt(sum((model - observation)^2, na.rm=na.rm)/(min(sum(!is.na(model)),sum(!is.na(observation)))))
}
