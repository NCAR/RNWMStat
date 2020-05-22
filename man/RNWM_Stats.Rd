\name{RNWM_Stats}
\alias{RNWM_Stats}
\title{RNWM_Stats (Estimate all the model evaluation statistics)}
\usage{
RNWM_Stats(dataFrame,bf_method='Eckhardt',qlimit=0,cor_method='pearson',cor_use='everything',na.rm=T)
}
\description{
Master function to estimate all the model evaluation statistics in 'RNWMStat' package. User should supply a dataframe of observation and model streamflow. Use column names 'mod' for model and 'obs' for observation.
}
\examples{

dataFrame<-data.frame(mod,obs)

nwm_validation_stats<-RNWM_Stats(dataFrame,bf_method='Eckhardt',qlimit=0,cor_method='pearson',cor_use='everything',na.rm=T)


}
\arguments{
\item{dataFrame}{A dataframe of observation and model streamflow. Use column names 'mod' for model and 'obs' for observation}

\item{bf_method}{Baseflow seperation methods (e.g., "Duncan" or "Lyne-Nathan" or "Eckhardt")}
\item{qlimit}{Minimum quickflow threshold. The StormEvents() function only consider events higher than this threshold while counting events}

\item{cor_method}{a character string indicating which correlation coefficient (or covariance) is to be computed. One of "pearson" (default), "kendall", or "spearman": can be abbreviated.}

\item{cor_use}{an optional character string giving a method for computing covariances in the presence of missing values. This must be (an abbreviation of) one of the strings "everything", "all.obs", "complete.obs", "na.or.complete", or "pairwise.complete.obs".}


}
\value{
Retun a list of all RNWMStat statistics.
}


