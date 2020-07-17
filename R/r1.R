# May 27, 2020; E. Towler
# Use Lamontagne github Rcode: https://github.com/JRLamontagne/Efficiency/blob/master/functions3BT.R
# to create a deriviate R code to calculate r1, Stedinger's (1981) lognormal estimator.
# These are calculated using the same input data as LBEm and LBE'm
# (whereby if there are zeros, an epsilon is added)
# Line to be added to /glade/u/home/arezoo/wrf_hydro/PyWrfHydroCalib/core/calib_workflow.R
# statCorR1 = r1(chrt.obj.nona.nozeros$mod, chrt.obj.nona.nozeros$obs)
# June12, 2020, added modification for Tau to be zero if (min(x) - tau <= 0) to avoid log(0) (happens with sites with many zeros)

r1 <- function(mod, obs){
   S <- mod#data[,2]   # Specify similuations
   O <- obs#data[,1]   # Specify observations
   # Get sample size of station data (how many observed and PRMS simulated streamflow data)
   n <- length(obs)#nrow(data)

   # Calculate correlation coefficient (calculate both r and r*=r1)
   r <- cor(S,O)   # Calculate correlation coefficient (rho) between simulations and observations


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
   if (tau_O<0 | tau_S<0){
     tau_O <- 0
     tau_S <- 0
   }

   # New added by Erin June 12, 2020, to deal with sites with many zeros, where tau ends up being min(O) or min(S); avoids log 0.
   if ( (min(O) - tau_O <= 0) | (min(S) - tau_S <= 0) ) {
     tau_O = 0 # As above, this means fitting LN2 at these sites.
     tau_S = 0 # As above, this means fitting LN2 at these sites.
   }

   u <- log(O-tau_O)      # Compute u based on observations of LN2 to use for LBE equation
   v <- log(S-tau_S)      # Compute v based on observations of LN2 to use for LBE equation
   # Estimator of correlation coefficient, r* (r1),  using Stedinger (1981) eqn 2
   # Stedinger (1981) eqn 3: estimator of variance of log of observations and simulations (CHECK!!!)
   s2_yOyS <- 1/n*sum((u-mean(u))*(v-mean(v)))
   s2_yO <- 1/n*sum((u-mean(u))^2)   # Stedinger (1981) eqn 3: estimator of variance of log of observations
   s2_yS <- 1/n*sum((v-mean(v))^2)    # Stedinger (1981) eqn 3: estimator of variance of log of simulations
   r1 <- (exp(s2_yOyS)-1)/sqrt((exp(s2_yO)-1)*(exp(s2_yS)-1)) # Stedinger (1981) eqn 2: estimator of correlation coefficient, r* (r1)
   return(r1)
}

