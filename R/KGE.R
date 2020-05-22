KGE <- function (model, observation, na.rm=TRUE, s.r=1, s.alpha=1, s.beta=1) {
  use <- if(na.rm) 'pairwise.complete.obs' else 'everything'
  r     <- cor(model, observation, use=use)
  alpha <- sd(model, na.rm=na.rm) / sd(observation, na.rm=na.rm)
  beta  <- mean(model, na.rm=na.rm) / mean(observation, na.rm=na.rm)
  eds = sqrt( (s.r*(r-1))^2 + (s.alpha*(alpha-1))^2 + (s.beta*(beta-1))^2 )
  kges = 1-eds
  return(kges)
}
