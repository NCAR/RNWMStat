\name{NRMSE}
\alias{NRMSE}
\title{NRMSE (Normalized Root Mean Squared Error)}
\usage{
NRMSE(model,observation,na.rm=TRUE)
}
\description{
Normalized Root Mean Squared Error- between model and observation
}
\examples{
data("USGS_Streamflow")

NRMSE(model,observation, na.rm=TRUE)

}
\arguments{
\item{model}{model streamflow time series}

\item{observation}{observation streamflow time series }

}
\value{
Retun normalized root mean squred error between model and observation.
}

