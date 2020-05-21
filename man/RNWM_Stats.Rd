\name{RNWM_Stats}
\alias{RNWM_Stats}
\title{RNWM_Stats (Estimate all the model evaluation statistics)}
\usage{
RNWM_Stats(dataFrame,bf_method='Eckhardt',qlimit=0,na.rm=T)
}
\description{
Master function to estimate all the model evaluation statistics in 'RNWMStat' package. User should supply a dataframe of observation and model streamflow. Use column names 'mod' for model and 'obs' for observation.
}
\examples{

dataFrame<-data.frame(mod,obs)

nwm_validation_stats<-RNWM_Stats(dataFrame,bf_method='Eckhardt',qlimit=0,na.rm=T)


}
\arguments{
\item{dataFrame}{A dataframe of observation and model streamflow. Use column names 'mod' for model and 'obs' for observation}

\item{bf_method}{Baseflow seperation methods (e.g., "Duncan" or "Lyne-Nathan" or "Eckhardt")}
\item{qlimit}{Minimum quickflow threshold. The StormEvents() function only consider events higher than this threshold while counting events}

}
\value{
Retun a list of all RNWMStat statistics.
}


