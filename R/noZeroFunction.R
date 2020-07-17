# Function that adds an epsilon to entire distribution if there are zeros
# Based on hydroGOF library, for NSE(... FUN=log, epsilon="Pushpalatha2012") based on:
# Pushpalatha, R., Perrin, C., Le Moine, N. and Andreassian, V. (2012).
# A review of efficiency criteria suitable for evaluating low-flow simulations.
# Journal of Hydrology, 420, 171-182. DOI: 10.1016/j.jhydrol.2011.11.055
# If there are no zeros in the distribution, returns original dataframe.
#noZeroFunction = function(mod, obs, period){
noZeroFunction = function(mod, obs){ #YL, removed period since not used
  zmin = min(mod, obs, na.rm = T) # Xia added na.rm=T on 2020/06/17
  if (zmin ==0) {
    #Following Push2012, though I don't get identical results, they are very close
    epsilon = mean(obs, na.rm=T)/100
    obs = obs + epsilon
    mod = mod + epsilon
  } # end if (zmin =0)
  #df = data.table::data.table(q_cms = mod, obs = obs, period = period)  
  df = data.table::data.table(q_cms = mod, obs = obs)  
  return(df)
} # end function

