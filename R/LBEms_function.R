# May 27, 2020; E. Towler
# Use Lamontagne github Rcode: https://github.com/JRLamontagne/Efficiency/blob/master/functions3BT.R
# to create a deriviate R code to calculate BNL3 parameters.
# This is called by the LBEms_function
# June12, 2020, added modification for Tau to be zero if (min(x) - tau <= 0) to avoid log(0) (happens with sites with many zeros)
#               and a jitter to the obs and mod if ther variance for a particular month is zero. (i.e., all values are zero or same value in a month)

parms <- function(data){
  S <- data[,2]   # Specify similuations
  O <- data[,1]   # Specify observations

  # New added by Erin June 12, 2020. For some sites, an entire month may be zeros
  # (or same value for entire month), which leads to var=0. Add small jitter.
  # Sometimes the first time run, a warning is issued, but does not seem to affect anything.
  # Warning message:
  # In stats::runif(length(x), -amount, amount) :
  #  '.Random.seed[1]' is not a valid integer, so ignored
  if (var(O)==0 | var(S)==0) {
    O = jitter(O)
    S = jitter(S)
  }

  # Get sample size of station data (how many observed and PRMS simulated streamflow data)
  n <- nrow(data)

  mu_O <- mean(O)  # Calculate mean of observations,mu_O (mean of column 1)
  ############## LBE ##############
  # Eq. 16
  # Calculate tau for observations and simulations
  if ((min(O)+max(O)-2*median(O))>0){
    tau_O <- (min(O)*max(O)-median(O)^2)/(min(O)+max(O)-2*median(O))     # Compute tau_O for each MC trial
  } else {
    tau_O <- 0   # Fit LN2 model instead of LN3 because Stedinger (1980) lower bound estimator cannot be computed
  }

  if ((min(S)+max(S)-2*median(S))>0){
    tau_S <- (min(S)*max(S)-median(S)^2)/(min(S)+max(S)-2*median(S))     # Compute tau_S for each MC trial
  } else {
    tau_S <- 0   # Fit LN2 model instead of LN3 because Stedinger (1980) lower bound estimator cannot be computed
  }

  #### ADDED 11.23.2019 to be consistent with Vogel's suggestion (11/1/2019) if tau values are negative set tau to zero. This means fitting LN2 at these sites.
  # If can't reproduce past results comment this out
 if (tau_O<0 | tau_S<0){
    tau_O <- 0
    tau_S <- 0
  }

  # New added by Erin June 12, 2020, to deal with sites with many zeros, where tau ends up being min(O) or min(S); avoids log 0.
  #if ( (min(O) - tau_O <= 0) | (min(S) - tau_S <= 0) ) { # commented out by Xia
  if ( (min(O) - tau_O <= 0.0000000001) | (min(S) - tau_S <= 0.0000000001) ) { # Erin: New June 17, 2020, when there are many zeros, sometimes this is just barely positive, but TauO should be zero in that case
    tau_O = 0 # As above, this means fitting LN2 at these sites.
    tau_S = 0 # As above, this means fitting LN2 at these sites.
  }

  u <- log(O-tau_O)      # Compute u based on observations of LN2 to use for LBE equation
  mu_u <- mean(u)
  var_u <- (1/(n-1))*sum((u-mean(u))^2)   # Population variance of
  sd_u <- sqrt(var_u)
  v <- log(S-tau_S)      # Compute v based on observations of LN2 to use for LBE equation
  mu_v <- mean(v)
  var_v <- (1/(n-1))*sum((v-mean(v))^2)   # Population variance of v
  sd_v <- sqrt(var_v)

  rho_log <- cor(u,v)

  # Estimator of correlation coefficient, r* (r1),  using Stedinger (1981) eqn 2
  # Stedinger (1981) eqn 3: estimator of variance of log of observations and simulations (CHECK!!!)
  s2_yOyS <- 1/n*sum((u-mean(u))*(v-mean(v)))
  s2_yO <- 1/n*sum((u-mean(u))^2)   # Stedinger (1981) eqn 3: estimator of variance of log of observations
  s2_yS <- 1/n*sum((v-mean(v))^2)    # Stedinger (1981) eqn 3: estimator of variance of log of simulations
  r1 <- (exp(s2_yOyS)-1)/sqrt((exp(s2_yO)-1)*(exp(s2_yS)-1)) # Stedinger (1981) eqn 2: estimator of correlation coefficient, r* (r1)

  # LN3 estimators
  Co_LBE <- (sqrt(exp(2*mean(u)+var_u)*(exp(var_u)-1)))/(tau_O+exp(mean(u)+var_u/2)) # LBE estimate of coefficient of variation of observations
  theta_LBE <- sqrt((exp(2*mean(v)+var_v)*(exp(var_v)-1))/(exp(2*mean(u)+var_u)*(exp(var_u)-1))) # LBE estimate of theta
  delta_LBE <- 1-(tau_S+exp(mean(v)+var_v/2))/(tau_O+exp(mean(u)+var_u/2)) # LBE estimate of delta

  LBE <- 2*theta_LBE*r1-theta_LBE^2-delta_LBE^2/Co_LBE^2 # Barber et al. 2019 eq. 16
  LBEprime <- 1-sqrt((-1*delta_LBE)^2+(theta_LBE-1)^2+(r1-1)^2)

  LN3parms <- cbind(mu_O, r1, Co_LBE, delta_LBE, theta_LBE, tau_O, tau_S, mu_u, mu_v, sd_u, sd_v, rho_log)

  return(LN3parms)
}



# May 19, 2020; E. Towler
# Use Lamontagne github Rcode  https://github.com/JRLamontagne/Efficiency/blob/master/LBEm.R
# to create a deriviate R code to calculate âLBEmâ and "LBE'm",
# From Lamontagne, J. R., Barber C, Vogel RM (in review).
# Improved Estimators of Model Performance Efficiency for Skewed Hydrologic Data,
# Water Resources Research.
# June 12, 2020. Generalized the code for cases where there is not any data in a month (e.g., frozen rivers in Jan/Feb.)

#LBEms_function <- function(mod, obs, period, calcDailyStats) { # Function reads in paired model and obs without any NAs.
LBEms_function <- function(mod, obs, period, monsample) { # Function reads in paired model and obs without any NAs.
# YL, 07-17-2020: replaced calcDailyStats with monsample in the argument list

  # Remove NA
  id<-!is.na(mod) &!is.na(obs); mod <- mod[id]; obs <- obs[id]; period<-period[id] 
  allData = data.frame("O"=obs, "S"= mod, "month" = period)

  mu_O_month <- matrix(nrow=1, ncol=12)
  Co_month <- matrix(nrow=1, ncol=12)
  delta_month <- matrix(nrow=1, ncol=12)
  theta_month <- matrix(nrow=1, ncol=12)
  rho_month <- matrix(nrow=1, ncol=12)
  tau_O_month<- matrix(nrow=1, ncol=12)
  tau_S_month<- matrix(nrow=1, ncol=12)
  mu_u_month<- matrix(nrow=1, ncol=12)
  mu_v_month<- matrix(nrow=1, ncol=12)
  sd_u_month<- matrix(nrow=1, ncol=12)
  sd_v_month<- matrix(nrow=1, ncol=12)
  rho_log_month<- matrix(nrow=1, ncol=12)
  # mixture moment estimators
  mu_mix_O_month<- matrix(nrow=1, ncol=12)
  var_mix_O_month<- matrix(nrow=1, ncol=12)
  mu_mix_S_month<- matrix(nrow=1, ncol=12)
  var_mix_S_month<- matrix(nrow=1, ncol=12)
  mu_mix_SO_month<- matrix(nrow=1, ncol=12)
  mu_u_mix<- matrix(nrow=1, ncol=12)
  mu_v_mix<- matrix(nrow=1, ncol=12)
  tau_O_mix<- matrix(nrow=1, ncol=12)
  tau_S_mix<- matrix(nrow=1, ncol=12)
  sd_u_mix<- matrix(nrow=1, ncol=12)
  sd_v_mix<- matrix(nrow=1, ncol=12)
  mu_mix_O_month_MC<- matrix(nrow=1, ncol=12)
  var_mix_O_month_MC_PART<- matrix(nrow=1, ncol=12)
  mu_mix_S_month_MC<- matrix(nrow=1, ncol=12)
  var_mix_S_month_MC_PART<- matrix(nrow=1, ncol=12)
  mu_mix_SO_month_MC<- matrix(nrow=1, ncol=12)
  rho_mix<- matrix(nrow=1, ncol=12)

  mu_mix_O_MC <- matrix(nrow=1, ncol=1)
  var_mix_O_MC <- matrix(nrow=1, ncol=1)
  mu_mix_S_MC <- matrix(nrow=1, ncol=1)
  var_mix_S_MC <- matrix(nrow=1, ncol=1)
  mu_mix_SO_MC <- matrix(nrow=1, ncol=1)
  theta_mix_MC <- matrix(nrow=1, ncol=1)
  delta_mix_MC <- matrix(nrow=1, ncol=1)
  Co_mix_MC <- matrix(nrow=1, ncol=1)
  r_mix_MC <- matrix(nrow=1, ncol=1)
  LBE_mix_MC <- matrix(nrow=1, ncol=1)
  LBEprime_mix_MC <- matrix(nrow=1, ncol=1)
  mu_mix_O <- matrix(nrow=1, ncol=1)
  var_mix_O <- matrix(nrow=1, ncol=1)
  mu_mix_S <- matrix(nrow=1, ncol=1)
  var_mix_S <- matrix(nrow=1, ncol=1)

  mu_mix_SO <- matrix(nrow=1, ncol=1)

  theta_mix <- matrix(nrow=1, ncol=1)
  delta_mix <- matrix(nrow=1, ncol=1)
  Co_mix <- matrix(nrow=1, ncol=1)
  Cs_mix <- matrix(nrow=1, ncol=1)
  r_mix <- matrix(nrow=1, ncol=1)
  LBE_mix <- matrix(nrow=1, ncol=1)
  LBEprime_mix <- matrix(nrow=1, ncol=1)
  rho <- matrix(nrow=1, ncol=1)

  #if (calcDailyStats) monsample <- 30  else monsample <- 10 
  for (m in 1:12){ # For calculations based on monthly data
    oneSite_month <- subset(allData, allData$month== m) # Pull data for individual site's month
    if (nrow(oneSite_month) > monsample ) { # added this line - Xia pointed out case where there's missing data for one month
      LN3params_month <- parms(oneSite_month[,1:2])  # Xia: 3-yr sample limit for 0 as we discussed? on 20200616
      mu_O_month[m] <- LN3params_month[1] # real space mean
      rho_month[m] <- LN3params_month[2]
      Co_month[m] <- LN3params_month[3]
      delta_month[m] <- LN3params_month[4]
      theta_month[m] <- LN3params_month[5]
      tau_O_month[m] <- LN3params_month[6]
      tau_S_month[m] <- LN3params_month[7]
      mu_u_month[m] <- LN3params_month[8]
      mu_v_month[m] <- LN3params_month[9]
      sd_u_month[m] <- LN3params_month[10]
      sd_v_month[m] <- LN3params_month[11]
      rho_log_month[m] <- LN3params_month[12]

      # Per Vogel's suggestion (11/1/2019) if tau values are negative set tau to zero. This means fitting LN2 at these sites.
      if (tau_O_month[m]<0 | tau_S_month[m]<0){
        tau_O_month[m] <- 0
        tau_S_month[m] <- 0
        LN2count_month <- LN2count_month + 1
      }

      # mixture moment estimators
      mu_mix_O_month[m] <- tau_O_month[m]+exp(mu_u_month[m]+sd_u_month[m]^2/2)
      var_mix_O_month[m] <- (exp(2*mu_u_month[m]+sd_u_month[m]^2)*(exp(sd_u_month[m]^2)-1))
      mu_mix_S_month[m] <- tau_S_month[m]+exp(mu_v_month[m]+sd_v_month[m]^2/2)
      var_mix_S_month[m] <- (exp(2*mu_v_month[m]+sd_v_month[m]^2)*(exp(sd_v_month[m]^2)-1))
      mu_mix_SO_month[m] <- (mu_mix_S_month[m]*mu_mix_O_month[m]+rho_month[m]*sqrt(var_mix_S_month[m])*sqrt(var_mix_O_month[m]))

      var_mix_S_month[m] <- (exp(2*mu_v_month[m]+sd_v_month[m]^2)*(exp(sd_v_month[m]^2)-1))
      mu_mix_SO_month[m] <- (mu_mix_S_month[m]*mu_mix_O_month[m]+rho_month[m]*sqrt(var_mix_S_month[m])*sqrt(var_mix_O_month[m]))
    } # end if (nrow(oneSite_month) > 0 ) { # added this line.
  } # end for (m in 1:12){ # For calculations based on monthly data

  # mixture moments from RAW data  - June 12, edited these to be mean( , na.rm =T), instead of 1/12*sum
  mu_mix_O  <- mean(mu_mix_O_month, na.rm =T)
  var_mix_O <- mean((var_mix_O_month+mu_mix_O_month^2), na.rm=T)-mu_mix_O^2
  mu_mix_S  <- mean(mu_mix_S_month, na.rm=T)
  var_mix_S <- mean((var_mix_S_month+mu_mix_S_month^2),na.rm=T)-mu_mix_S^2

  mu_mix_SO <- mean(mu_mix_SO_month, na.rm=T)

  theta_mix <- sqrt(var_mix_S)/sqrt(var_mix_O)
  delta_mix <- 1-mu_mix_S/mu_mix_O
  Co_mix <- sqrt(var_mix_O)/mu_mix_O
  Cs_mix <- sqrt(var_mix_S)/mu_mix_S
  #r1_mix[i] <- (1/nrow(oneSite)*sum(oneSite[,3]*oneSite[,4])-mu_mix_O[i]*mu_mix_S[i])/(sqrt(var_mix_O[i]*var_mix_S[i]))
  r_mix <- (mu_mix_SO-mu_mix_O*mu_mix_S)/(sqrt(var_mix_O*var_mix_S))

  LBE_mix <- 2*theta_mix*r_mix-theta_mix^2-delta_mix^2/Co_mix^2
  LBEprime_mix <- 1-sqrt(delta_mix^2+(theta_mix-1)^2+(r_mix-1)^2)
  LBEms = c(LBE_mix, LBEprime_mix)
  return(LBEms)
} #end
                                                      

