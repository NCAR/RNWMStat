\name{RMSE}
\alias{RMSE}
\title{RMSE (Root Mean Squared Error)}
\usage{
RMSE(model,observation,na.rm=TRUE)
}
\description{
Root Mean Squared Error- between model and observation
}
\examples{
data("USGS_Streamflow")

RMSE(mod,obs)

}
\arguments{
\item{model}{model streamflow time series}

\item{observation}{observation streamflow time series }

}
\value{
Retun root mean squred error between model and observation.
}


